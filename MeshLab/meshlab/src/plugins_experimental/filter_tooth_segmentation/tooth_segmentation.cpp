/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/


#include "tooth_segmentation.h"
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/create/platonic.h> // �������򼸺���״�Ļ���
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <vcg/complex/algorithms/update/geodesic.h>
#include <vcg/complex/algorithms/update/nring.h>
#include <vcg/math/disjoint_set.h>
#include <vcg/simplex/face/pos.h>
#include <wrap\io_trimesh\export_stl.h>
#include <QtScript>
#include <meshlab/glarea.h>
#include <iostream>
#include <fstream>

#include <stdio.h>
#define PI 3.14159265


using namespace std;
using namespace vcg;


ToothSegmentationPlugin::ToothSegmentationPlugin(void)
{
    typeList
        << FP_MESH_SEGMENTATION
        ;

    FilterIDType tt;

    foreach(tt, types())
        actionList << new QAction(filterName(tt), this);

}

ToothSegmentationPlugin::FilterClass ToothSegmentationPlugin::getClass(QAction * a)
{
	switch (ID(a))
	{
	case FP_MESH_SEGMENTATION:  return MeshFilterInterface::Smoothing;
	default: assert(0); return MeshFilterInterface::Generic;
	}
}

QString ToothSegmentationPlugin::filterName(FilterIDType filter) const
{
    switch (filter)
    {
	case FP_MESH_SEGMENTATION:  return QString("Automatic Tooth Segmentation");
	default: assert(0);
    }
    return tr("error!");
}

QString ToothSegmentationPlugin::filterInfo(FilterIDType filterID) const
{
    switch (filterID)
    {
	case FP_MESH_SEGMENTATION:  return QString("Automatic tooth segmentation on dental meshes using morphological operators.\n"
		" Warrning: This version works only with closed meshes.");
    default                                  : assert(0);
    }

    return QString();
}

// this function builds and intializes with the default values (that can depend on the current mesh or selection)
// the list of parameters that a filter requires.
// return
//		true if has some parameters
//		false is has no params
void ToothSegmentationPlugin::initParameterSet(QAction * action, MeshModel & m, RichParameterSet & parlst)
{   
    QStringList methods;
    QStringList loopWeightLst;

    switch(ID(action))
    {
	case FP_MESH_SEGMENTATION:
		parlst.addParam(new RichFloat("ThresholdValue",
			0.45f,
			"Threshold of Curvature",
			"Initial feature region according mean curvature."));
		parlst.addParam(new RichFloat("MaxDist",
			3.5f,
			"Max Neighborhood Dist",
			"Find middle point according Max Neighborhood Dist."));
		break;

    default:
        break;
    }
}

bool ToothSegmentationPlugin::applyFilter(QAction * filter, MeshDocument & md, RichParameterSet & par, vcg::CallBackPos * cb)
{
	MeshModel &m = *md.mm();

	switch (ID(filter))
	{
	case FP_MESH_SEGMENTATION:
	{
		MeshModel &m = *(md.mm());

		//Enable curvature and adjacency relations needed to compute a curvature.
		m.updateDataMask(MeshModel::MM_VERTCURV | MeshModel::MM_VERTCURVDIR);
		m.updateDataMask(MeshModel::MM_VERTCOLOR | MeshModel::MM_FACECOLOR);		
		m.updateDataMask(MeshModel::MM_VERTMARK | MeshModel::MM_FACEMARK);
		m.updateDataMask(MeshModel::MM_VERTFACETOPO | MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTQUALITY);

		// IsB() ��ǵĸ���
		vcg::tri::UpdateFlags<CMeshO>::FaceBorderFromVF(m.cm);
				
		setVFB = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<char>(m.cm, string("setFeatureBorderVertex"));
		marked = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(m.cm, string("marked"));
		curvature = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<float>(m.cm, string("curvature"));
		disc = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(m.cm, string("disc"));
		center = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(m.cm, string("center"));
		FSetNum = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(m.cm, string("FeatureSetNumber"));
		disth = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<float>(m.cm, std::string("DistFromFeat"));
		CpointOwner = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<int>(m.cm, std::string("CpointOwner"));
		divfval = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<float>(m.cm, std::string("DivergenceFieldValue"));
		source = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<CMeshO::VertexPointer>(m.cm, std::string("source"));
		source1 = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<CMeshO::VertexPointer>(m.cm, std::string("source1"));

		vcg::tri::Clean<CMeshO>::RemoveUnreferencedVertex(m.cm);
		vcg::tri::Allocator<CMeshO>::CompactVertexVector(m.cm);
		
		vcg::tri::UpdateNormal<CMeshO>::NormalizePerVertex(m.cm);
		vcg::tri::UpdateCurvatureFitting<CMeshO>::computeCurvature(m.cm);	//k1k2���ʹ���				
		vcg::tri::UpdateQuality<CMeshO>::VertexFromMeanCurvatureDir(m.cm);  //ƽ�����ʼ��㲢��������
				
		
		// ����������ֵɸѡ����
		if ( CurvGroup(m, par) == false)
			return false;

		// ��ʼ��������ʶ
		CMeshO::VertexIterator vi;
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) 
		{	//all vertices start as set V for Vertex and unmarked
			setVFB[vi] = 'V';
			
			marked[vi] = 0;
			disc[vi] = 0;		// ����������, �����ķ� complex �㶼�� disk
			center[vi] = 0;		// �������� center vertex �ı��, ��Χ���� feature vertex ʱ���Ǹõ���� center, ���Ϊ 1
			
			//all vertices start as -1 FSetNum
			FSetNum[vi] = -1;
			//none have been visited to connect C points, so set CpointOwner to -1
			CpointOwner[vi] = -1;

			curvature[vi] = (*vi).Q();
		}
		for (vector<CMeshO::VertexPointer>::iterator iterFP = FeaturePoints.begin();
			iterFP != FeaturePoints.end(); ++iterFP)
		{
			// feature region, less than defined min curvature
 			setVFB[(*iterFP)] = 'F';
		}		
		
		FilterNoise(m, 0.01);
				
		// morphologic operation
		meshDilate(m, 2);		
		meshErode(m, 1);		
		meshSkeletonize(m);		
		meshPrune(m, 4);  // ���� 4 �μ�������

		FilterNoise(m, 0.02);		
		
		vector<CVertexO*> ConnectPoints = FindConnectPoints(m);
		vector<CMeshO::VertexPointer> MiddlePoints = FindMidPoints(m, par, ConnectPoints);
		//CalculateDistanceField(m, ConnectPoints);
		AddPathtoConnect(m, MiddlePoints);
		showSegments1st(m);
		RemoveLines(m);

		//SaveAsSTL(m);
		
		ShowFeature(m);

		/*
		for (vector<CVertexO*>::iterator iterV  = ConnectPoints.begin(); iterV != ConnectPoints.end(); ++iterV)
		{// ��ʾ���ӵ�
		(*iterV)->SetS();
		}
		for (vector<CMeshO::VertexPointer>::iterator iterV  = MiddlePoints.begin();
		iterV != MiddlePoints.end(); ++iterV)
		{// ��ʾ�м��
		(*iterV)->C() = vcg::Color4b(vcg::Color4b::Blue);
		}*/

		
		// delete attibute by name
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<int>(m.cm, marked);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<float>(m.cm, disth);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<char>(m.cm, setVFB);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<float>(m.cm, divfval);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<int>(m.cm, FSetNum);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<int>(m.cm, CpointOwner);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<CMeshO::VertexPointer>(m.cm, source);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<float>(m.cm, curvature);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<int>(m.cm, disc);
		vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<int>(m.cm, center);
		return true;
	} break;
	}
	return true;
}

int ToothSegmentationPlugin::postCondition(QAction * filter) const
{
    switch (ID(filter))
    {
    case FP_MESH_SEGMENTATION:
		return MeshModel::MM_VERTCOLOR | MeshModel::MM_VERTQUALITY | MeshModel::MM_VERTNORMAL;
    default                  : return MeshModel::MM_UNKNOWN;
    }
}

void ToothSegmentationPlugin::SaveAsSTL(MeshModel &m)
{
	CMeshO::VertexIterator vi;
	// �����ܹ��ж�������ɫ
	std::vector<vcg::Color4b> colorType;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		vcg::Color4b colorFliter = (*vi).C();
		vcg::Color4b colorTemp;
		if (!colorType.empty())
		{
			std::vector<vcg::Color4b>::iterator iter;
			for (iter = colorType.begin(); iter != colorType.end(); ++iter)
			{
				if (*iter == colorFliter) break;
			}
			if (iter == colorType.end()) colorType.push_back(colorFliter);

		}
		else
		{
			colorType.push_back(colorFliter);
		}
	}

	// ��ÿ���������ڵ���, ɾ����
	vcg::Color4b colorV = vcg::Color4b::LightGray;
	// if (!colorType.empty()) colorV = colorType[4];
	std::vector<CVertexO*> vecFaceV;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		vcg::Color4b colorFliter = (*vi).C();

		if (!vi->IsD() && (colorFliter == colorV || setVFB[vi] == 'F'))
		{
			vcg::face::VFIterator<CFaceO> vfi(&(*vi)); //initialize the iterator tohe first face
			for (; !vfi.End(); ++vfi)
			{
				CFaceO* f = vfi.F();
				// ...do something with face f
				if (!f->IsD())
				{
					//f->V(0)->SetD();
					//f->V(1)->SetD();
					//f->V(2)->SetD();
					vecFaceV.push_back(f->V(0));
					vecFaceV.push_back(f->V(1));
					vecFaceV.push_back(f->V(2));
					f->SetD();
				}
			}
		}
	}
	for (auto iter = vecFaceV.begin(); iter != vecFaceV.end(); ++iter) if (!(*iter)->IsD())
	{
		(*iter)->SetD();
	}

	vcg::tri::Allocator<CMeshO>::CompactFaceVector(m.cm);
	vcg::tri::Allocator<CMeshO>::CompactVertexVector(m.cm);
	vcg::tri::Allocator<CMeshO>::CompactEveryVector(m.cm);
	
	vcg::tri::io::ExporterSTL<CMeshO>::Save(m.cm, "delete_teeth_1.stl");
}

bool ToothSegmentationPlugin::CurvGroup(MeshModel &m, RichParameterSet &par)
{
	if (vcg::tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m.cm) >0 || vcg::tri::Clean<CMeshO>::CountNonManifoldVertexFF(m.cm) >0)
	{
		this->errorMessage = "Mesh has some not 2-manifold faces, cannot compute principal curvature directions"; // text
		return false; // can't continue, mesh can't be processed
	}

	CMeshO::VertexIterator vi;	

	FeaturePoints.clear();
	
	CurvThreshold = par.getFloat("ThresholdValue");
	/*FILE *fp;*/
	/*fp = fopen("test.txt", "w");*/
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		if (!(*vi).IsD())
		{	//float fTemp = 2.0 / (1 + exp(-5.0 * (*vi).Q())) - 1.0;
			float fTemp = (*vi).Q();
			/*fprintf(fp, "%4.8f\n", fTemp);*/
			if (fTemp < -1.0 * CurvThreshold) FeaturePoints.push_back(&*vi);			
		}
	}
	/*fclose(fp);*/
	return true;	
}

void ToothSegmentationPlugin::ShowFeature(MeshModel &m)
{
	vcg::tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
		
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (setVFB[vi] == 'F')
		{
			(*vi).C() = vcg::Color4b::Red;			
		}
	}
}

void ToothSegmentationPlugin::meshDilate(MeshModel &m, int n_ring)
{
	CMeshO::VertexIterator vi;
	//dilate 1-ring neighbor
	for (int k = 0; k < n_ring; k++) {
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (setVFB[vi] == 'V') {
				(*vi).SetV();				
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
				OneRing.insertAndFlag1Ring(&(*vi)); // ���������� vi ���ڵĵ�	
				(*vi).ClearV();
				for (int i = 0; i < OneRing.allV.size(); i++) {
					CVertexO* tempVertexPointer = OneRing.allV.at(i);
					//if the neighbor of the current vertex is not a feature, mark the neighbor to be one
					//if it is negative curvature
					if (setVFB[&(*tempVertexPointer)] == 'F') {
						//mark it to be 'N'						
						marked[vi] = 1; break;
					}
				}				
			}
		}
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			//if it was marked, change it
			if (marked[vi] == 1) {
				setVFB[vi] = 'F';
				marked[vi] = 0;
			}
		}
	}	
}

void ToothSegmentationPlugin::meshErode(MeshModel &m, int n_ring)
{
	CMeshO::VertexIterator vi;
	//erode 1-ring neighbor
	for (int k = 0; k < n_ring; k++) {
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (setVFB[vi] == 'F') {
				(*vi).SetV();
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
				OneRing.insertAndFlag1Ring(&(*vi));		
				(*vi).ClearV();			
				for (int i = 0; i < OneRing.allV.size(); i++) {
					CVertexO* tempVertexPointer = OneRing.allV.at(i);
					//if the neighbor of the current vertex is not a feature, we will erode this vertex
					if (setVFB[&(*tempVertexPointer)] == 'V') {
						marked[vi] = 1; break;
					}
				}				
			}
		}
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			//if it was marked to not erode, keep it
			if (marked[vi] == 1) {
				setVFB[vi] = 'V';
				marked[vi] = 0;
			}
		}
	}
}

void ToothSegmentationPlugin::meshSkeletonize(MeshModel &m)
{
	CMeshO::VertexIterator vi;
	int debugCounter = 1;
	bool changes = false;
	//�������ڱ䶯�ͼ���ִ�� skeletonize
	do {
		changes = false;
		//����������, �ж��� centers or disks
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (setVFB[vi] == 'F') {
				//check for centers and disks
				(*vi).SetV();
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
				OneRing.insertAndFlag1Ring(&(*vi));
				(*vi).ClearV();
				bool neighborsAllF = true;
				for (int i = 0; i < OneRing.allV.size(); i++) {
					CVertexO* tempVertexPointer = OneRing.allV.at(i);
					//vi����㲻��������, vi�Ͳ��� center
					if (setVFB[&(*tempVertexPointer)] != 'F') {
						neighborsAllF = false;
					}
				}
				if (neighborsAllF) {
					//vi���Ϊ center
					center[vi] = 1;
					//vi�������Ϊ disc
					for (int i = 0; i < OneRing.allV.size(); i++) {
						CVertexO* tempVertexPointer = OneRing.allV.at(i);
						disc[&(*tempVertexPointer)] = 1;
					}
				}
			}//in set F
		}//go through vertices
		 //������center ��Ϊdisc�Ķ���, ����complexity, ɾ����complex
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (center[vi] == 0 && disc[vi] == 1) {
				//check complexity
				(*vi).SetV();
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
				OneRing.insertAndFlag1Ring(&(*vi));	
				(*vi).ClearV();		
				int numberOfChanges = 0;    // ��¼�õ�� complexity
				char iVertex;
				char iPlus1Vertex;
				CVertexO* tempVertexPointer;
				for (int i = 0; i < OneRing.allV.size(); i++) {
					if (i < (OneRing.allV.size() - 1)) {
						tempVertexPointer = OneRing.allV.at(i);
						iVertex = (setVFB[&(*tempVertexPointer)]);
						tempVertexPointer = OneRing.allV.at(i + 1);
						iPlus1Vertex = (setVFB[&(*tempVertexPointer)]);
					}
					else {
						tempVertexPointer = OneRing.allV.at(i);
						iVertex = (setVFB[&(*tempVertexPointer)]);
						tempVertexPointer = OneRing.allV.at(0);
						iPlus1Vertex = (setVFB[&(*tempVertexPointer)]);
					}
					//change noted
					if (iVertex != iPlus1Vertex) {
						numberOfChanges++;
					}
				}
				//�� complex, ���������� F ��ɾ��
				if (numberOfChanges < 4) {
					//delete from feature
					setVFB[vi] = 'V';
					center[vi] = 0;
					disc[vi] = 0;
					changes = true;
				}
			}//not a center, yes a disk
		}//go through vertices

		 //************ һ�� skeletonize ����,ɾ�� disk, center ���
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (setVFB[vi] == 'F') {
				center[vi] = 0;
				disc[vi] = 0;
			}
		}

		++debugCounter;
	} while (changes == true);//loop if changes
}

void ToothSegmentationPlugin::meshPrune(MeshModel &m, int pruneLength)
{
	CMeshO::VertexIterator vi;
	//pruneLength �� prune ����
	for (int pruneIter = 0; pruneIter < pruneLength; pruneIter++) {
		//check complexity if the vertices is 'F'
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (setVFB[vi] == 'F') {
				//check complexity
				(*vi).SetV();
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
				OneRing.insertAndFlag1Ring(&(*vi));	
				(*vi).ClearV();
				int numberOfChanges = 0;
				char iVertex;
				char iPlus1Vertex;
				CVertexO* tempVertexPointer;
				for (int i = 0; i < OneRing.allV.size(); i++) {
					if (i < (OneRing.allV.size() - 1)) {
						tempVertexPointer = OneRing.allV.at(i);
						iVertex = (setVFB[&(*tempVertexPointer)]);
						tempVertexPointer = OneRing.allV.at(i + 1);
						iPlus1Vertex = (setVFB[&(*tempVertexPointer)]);
					}
					else {
						tempVertexPointer = OneRing.allV.at(i);
						iVertex = (setVFB[&(*tempVertexPointer)]);
						tempVertexPointer = OneRing.allV.at(0);
						iPlus1Vertex = (setVFB[&(*tempVertexPointer)]);
					}
					//change noted
					if (iVertex != iPlus1Vertex) {
						numberOfChanges++;
					}
				}
				//not complex, delete it from feature
				if (numberOfChanges < 4) {
					//prune from feature
					setVFB[vi] = 'V';					
				}
			}//checking current vertex
		}//go through vertices
	}//number of times we prune
}

void ToothSegmentationPlugin::FilterNoise(MeshModel &m, float deletePercentOfF)
{
	CMeshO::VertexIterator vi;

	// ***** ���� Disjoint Sets �ṹ���������㰴���ϻ���, ÿ��������������������һ�� set
	// ��ʼ������, ÿ�������㵥����Ϊһ�� set�� ֮���� Union ����
	vcg::DisjointSet<CVertexO>* ptrDset = new vcg::DisjointSet<CVertexO>();
	// Initialize the Disjoint sets of 'F' and delete any that are small
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (setVFB[vi] == 'F') {
			//put the vertex in the set D
			ptrDset->MakeSet(&(*vi));   // ��ʼʱ, ÿ�����㵥����Ϊһ�� set
		}
	}
	
	// ***** ���������������� Union ����, �鲢��ͬһ������
	int blueCounter = 0;    // ��¼ feature points ����, Ϊ�����ɾ���������ο�
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//if it is a feature, then find 1-ring and use Union to merge sets
		if (setVFB[vi] == 'F') {
			blueCounter++;			
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));	
			(*vi).ClearV();
			for (int i = 0; i < OneRing.allV.size(); i++) {
				CVertexO* tempVertexPointer = OneRing.allV.at(i);
				//if the neighbor of the current vertex is a feature then merge sets
				if (setVFB[&(*tempVertexPointer)] == 'F') {
					//if they are not already in the same set, merge them
					if (ptrDset->FindSet(&(*vi)) != ptrDset->FindSet(tempVertexPointer)) {
						ptrDset->Union(&(*vi), tempVertexPointer); // ����������ϲ�ͬ���ͺϲ�
					}
				}
			}
		}
	}
	
	// ***** ����洢��ͬ�ļ��ϵ� allSetsWithVertices ��
	// ���Ȳ��Ұ����Ϸ�������ݳ�ʼ��, �ٳ�ʼ�� allSetsWithVertices ��������, 
	// ���ͨ�� Find ��������ͬһ���ϵĵ�, ������

	// Find out how many distinct sets there are
	int setCounter = 0;
	vector <CVertexO*> parentVertexofSet;
	typedef vector < vector <CVertexO*> > matrix;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//if it is in set 'F', then find what set it belongs to
		if (setVFB[vi] == 'F') {
			if (ptrDset->FindSet(&(*vi)) == &(*vi)) {
				//parent of a new set
				parentVertexofSet.push_back(&(*vi));
				setCounter++;
			}
		}
	}
	matrix allSetsWithVertices(setCounter, vector<CVertexO*>(0));
	// Create a vector of vectors containing our sets / vertices
	vector<CVertexO*>::iterator itVect;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		// if it is in set 'F', then find what set it belongs to
		if (setVFB[vi] == 'F') {
			// find offset in parent vector--that is the vertex set number
			itVect = find(parentVertexofSet.begin(), parentVertexofSet.end(), ptrDset->FindSet(&(*vi)));
			if (itVect != parentVertexofSet.end()) {
				int offset = std::distance(parentVertexofSet.begin(), itVect);
				//store the set number in the vertex
				FSetNum[vi] = offset;
				allSetsWithVertices[offset].push_back(&(*vi));
				// allSetsWithVertices �а���洢��ÿ����������ĵ㼯
			}
		}
	}
	
	// ***** ɾ��С�������㼯
	delete ptrDset;
	//delete sets of D if the number of vertices is less than deletePercentOfF
	//���� deletePercentOfF = 0.01;
	for (int pos = 0; pos < allSetsWithVertices.size(); pos++)
	{
		if (allSetsWithVertices[pos].size() <= deletePercentOfF*blueCounter) {
			for (int i = 0; i < allSetsWithVertices[pos].size(); i++) {
				// (*(allSetsWithVertices[pos][i])).C() = vcg::Color4b::White;
				setVFB[(*(allSetsWithVertices[pos][i]))] = 'V';
			}
			allSetsWithVertices[pos].clear();
		}
	}	
}

void ToothSegmentationPlugin::showSegments1st(MeshModel &m)
{
	// ============
	//�Բ�ͬ�� segments ��ɫ, ��ɾ��������ɫ��ͬ�� feature lines 
	float deletePercentOfBiggestSegment = 0.25;
	CMeshO::VertexIterator vi;

	// ***** �Է� 'F' �㴴�� Disjoint Sets
	vcg::DisjointSet<CVertexO>* DisjointSetsToColor = new vcg::DisjointSet<CVertexO>();
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//Every vertex not in set 'F' starts as disjoint set
		if (setVFB[vi] != 'F' /*&& (!(*vi).IsB())*/) {
			//put the vertex in the set D
			DisjointSetsToColor->MakeSet(&(*vi));
		}
	}

	// ***** Merge neighboring Sets
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		// if it is not a feature, then find 1-ring and merge sets
		if (setVFB[vi] != 'F' /*&& (!(*vi).IsB())*/) {
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));	
			(*vi).ClearV();
			for (int i = 0; i < OneRing.allV.size(); i++) {
				CVertexO* tempVertexPointer = OneRing.allV.at(i);
				// if the neighbor of the current vertex is not a feature then merge sets
				if (setVFB[&(*tempVertexPointer)] != 'F') {
					// if they are not already in the same set, merge them
					if (DisjointSetsToColor->FindSet(&(*vi)) != DisjointSetsToColor->FindSet(tempVertexPointer)) {
						DisjointSetsToColor->Union(&(*vi), tempVertexPointer);
					}
				}
			}
		}
	}   // ��ÿ����� 1-ring �ڽӵ� Union ������, �����ÿ����������������ڵĵ���ͬһ����

		// ***** ���� segments �ĸ���
	int segmentSetCounter = 0;
	vector <CVertexO*> parentVertexofDisjointSets;
	// If a vertex is the parent of a DisjointSet, it is a new Segment
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		// if it is not in set 'F', then find what set it belongs to
		if (setVFB[vi] != 'F') {
			if (DisjointSetsToColor->FindSet(&(*vi)) == &(*vi)) {
				// parent of a new set
				parentVertexofDisjointSets.push_back(&(*vi));
				segmentSetCounter++;
			}
		}
	}

	// ===== find biggest set
	// ***** ����洢��ͬ segment �Ķ���
	typedef vector < vector <CVertexO*> > crmatrix;
	crmatrix vectorOfDisjointSetsWithVertices(segmentSetCounter, vector<CVertexO*>(0));
	vector<CVertexO*>::iterator itVect;
	// ************* �� vector ������ṹ�洢��ͬ��ļ��� sets / vertices
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//if is 'F', then find what set it belongs to
		if (setVFB[vi] != 'F') {
			//�� find �������� vi ������ set �� matrix �е�λ��
			itVect = find(parentVertexofDisjointSets.begin(), parentVertexofDisjointSets.end(), DisjointSetsToColor->FindSet(&(*vi)));
			if (itVect != parentVertexofDisjointSets.end()) {
				int offset = std::distance(parentVertexofDisjointSets.begin(), itVect);
				// vi ��¼����Ӧ�� set ��
				FSetNum[vi] = offset;
				vectorOfDisjointSetsWithVertices[offset].push_back(&(*vi));
			}
		}
	}

	//////////////////////////////////////////////////////////////
	// ***** For each set, sum the positive curvature and store in new vector
	vector<float> vectorOfSetPosCurvatureSums(vectorOfDisjointSetsWithVertices.size(), 0);
	for (int pos = 0; pos < vectorOfDisjointSetsWithVertices.size(); pos++)
	{
		float tempSumOfCurvature = 0.0;
		for (int numInSet = 0; numInSet < vectorOfDisjointSetsWithVertices[pos].size(); numInSet++) {
			CVertexO* tempVertex = vectorOfDisjointSetsWithVertices[pos][numInSet];
			float tempCurvature = curvature[tempVertex];
			// if it's positive, add it up
			if (tempCurvature > 0.0) {
				tempSumOfCurvature = tempSumOfCurvature + tempCurvature;
			}
		}
		vectorOfSetPosCurvatureSums[pos] = tempSumOfCurvature;
	}   // �ⲿ�ְ�ÿ�����򶥵������ֵΪ�����ܺʹ�����,

		// ***** Find the set with the most vertices
	int maxInSet = 0;
	int maxOffset = 0;
	int secondBiggest = 0;  // �ڶ��� segment �����ݵ�����
	int secondBiggestOffset = 0;
	for (int pos = 0; pos < vectorOfDisjointSetsWithVertices.size(); pos++)
	{
		if (vectorOfDisjointSetsWithVertices[pos].size() > maxInSet) {
			secondBiggest = maxInSet;
			secondBiggestOffset = maxOffset;
			maxInSet = vectorOfDisjointSetsWithVertices[pos].size();
			maxOffset = pos;
		}
	}
	int numBigSets = 0;
	int numLittleSets = 0;
	for (int pos = 0; pos < vectorOfDisjointSetsWithVertices.size(); pos++)
	{
		if (pos != maxOffset) {
			if (vectorOfDisjointSetsWithVertices[pos].size() >(deletePercentOfBiggestSegment * secondBiggest)) {
				numBigSets++;   // ������������ 25% ��Ϊ������
			}
			else {
				numLittleSets++;
			}
		}
	}

	// ***** color the different sets
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (setVFB[vi] == 'F') {
			(*vi).C() = vcg::Color4b::Red;
		}
		// if it is not in set 'F', then find what set it belongs to
		if (setVFB[vi] != 'F') {
			if (FSetNum[vi] == maxOffset) {
				(*vi).C() = vcg::Color4b::LightGray;    // �����������������δ�ָ��
			}
			else {
				// color ramp
				int ScatterSize = vectorOfDisjointSetsWithVertices.size();
				vcg::Color4b BaseColor = vcg::Color4b::Scatter(ScatterSize, FSetNum[vi], .3f, .9f);
				(*vi).C() = BaseColor;
			}
		}
	}	
}

void ToothSegmentationPlugin::CalculateDistanceField(MeshModel &m)
{
	CMeshO::VertexIterator vi;
	float distanceConstraint = 0.8f;    // ������ Feature Line 0.5mm ��Χ�ڵľ��볡


										// ***** Mark the borders after deleting features
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//if it is a feature, then find 1-ring and mark set 'B'
		if (setVFB[vi] == 'F') {
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));
			(*vi).ClearV();
			for (int i = 0; i < OneRing.allV.size(); i++) {
				CVertexO* tempVertexPointer = OneRing.allV.at(i);
				if (setVFB[&(*tempVertexPointer)] == 'V') {
					// if the neighbor of the 'F' vertex is just a vertex it is on the border
					setVFB[&(*tempVertexPointer)] = 'B';
				}
			}
		}
	}

	// ***** Estimate distance from feature region
	vcg::tri::Geodesic<CMeshO> g;
	//create a vector of starting vertices in set 'B' (border of feature)
	std::vector<CVertexO*> fro;
	bool ret;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
		if (setVFB[vi] == 'B') {
			fro.push_back(&(*vi));
		}

	float dist_thr = distanceConstraint;
	if (!fro.empty()) {
		// ���㵽 fro �����е������ dist_thr ��Χ�ڵĲ�ؾ���
		ret = g.DistanceFromFeature(m.cm, fro, dist_thr);
	}
	float maxdist = distanceConstraint;
	float mindist = 0;
	// the distance is now stored in the quality, transfer the quality to the distance
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if ((*vi).Q() < distanceConstraint) {
			if (setVFB[vi] == 'F') {
				disth[vi] = -((*vi).Q());
			}
			if (setVFB[vi] == 'B') {
				disth[vi] = 0;
			}
			if (setVFB[vi] == 'V') {
				disth[vi] = (*vi).Q();
				setVFB[vi] = 'B';
			}
		}
		else {
			disth[vi] = std::numeric_limits<float>::max();
		}
	}

	/////////////////////////////////////////////
	// ***** filter the distances
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		// new distance to vertex is average of neighborhood distances
		if (disth[vi] < maxdist) {
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));
			(*vi).ClearV();
			float numberOfNeigbors = 0.0;
			float sumOfDistances = 0.0;
			for (int i = 0; i < OneRing.allV.size(); i++) {
				CVertexO* tempVertexPointer = OneRing.allV.at(i);
				float tempDist = disth[(*tempVertexPointer)];
				if (tempDist != std::numeric_limits<float>::max()) {
					sumOfDistances = sumOfDistances + tempDist;
					numberOfNeigbors++;
				}
			}
			float inverseNumNeighbors = 1.0 / numberOfNeigbors;
			disth[vi] = inverseNumNeighbors * sumOfDistances;
		}
	}   //////////////////////////////////////////////////////
}

void ToothSegmentationPlugin::CalculateDistanceField(MeshModel &m, std::vector<CVertexO*> vectCPoints)
{
	for (vector<CVertexO*>::iterator iterV = vectCPoints.begin(); iterV != vectCPoints.end(); ++iterV)
	{// ��ʾ���ӵ�
	 //(*iterV)->SetS();
		setVFB[(*iterV)] = 'B';
	}
}

vector<CVertexO *> ToothSegmentationPlugin::FindConnectPoints(MeshModel &m)
{
	CMeshO::VertexIterator vi;
	float connectingPointAngle = 240.0f;
	/**************** Connecting points by angle *****************/
	//find connecting points and color them yellow
	vector <CVertexO*> vectCPoints;
	int CPointCounter = 0;
	int vertexCounter = 0;
	//Find Connecting Points
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (setVFB[vi] == 'F') {
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));
			(*vi).ClearV();
			char currentSet = 'a';
			char initialSet = 'a';
			int numberOfChanges = 0;
			float numberOfNeighbors = 0.0;
			float sumOfDivergence = 0.0;
			float angle = 0.0;
			bool addAngle = false;
			CVertexO* firstVertex;
			CVertexO* previousVertex;
			char previousSet = 'a';
			char firstSet = 'a';
			int borderNeighborCounter = 0;
			int featureNeighborCounter = 0;
			for (int i = 0; i < OneRing.allV.size(); i++) {
				CVertexO* tempVertexPointer = OneRing.allV.at(i);
				char tempSetOfThisNeighbor = (setVFB[&(*tempVertexPointer)]);
				if (tempSetOfThisNeighbor == 'F') {
					featureNeighborCounter++;
					//don't add angle, we are still in 'F'
					addAngle = false;
				}
				if (tempSetOfThisNeighbor == 'B') {
					borderNeighborCounter++;
					//add angle, we are still in 'F'
					addAngle = true;
				}
				if (i == 0) {
					//initialize the starting set
					currentSet = setVFB[&(*tempVertexPointer)];
					initialSet = setVFB[&(*tempVertexPointer)];
					firstVertex = tempVertexPointer;
					firstSet = setVFB[&(*tempVertexPointer)];
					previousVertex = tempVertexPointer;
					previousSet = setVFB[&(*tempVertexPointer)];
				}
				if (tempSetOfThisNeighbor != currentSet) {
					//changed sets, update the numberOfChanges
					currentSet = setVFB[&(*tempVertexPointer)];
					numberOfChanges++;
				}
				//check if last point is in same set as first point
				if (i == (OneRing.allV.size() - 1)) {
					if (tempSetOfThisNeighbor != initialSet) {
						numberOfChanges++;
					}
					if (firstSet == 'F' && currentSet == 'F') {
						//don't add the angle
					}
					else {
						//calculate and add the angle between tempVertexPointer and firstVertex
						vcg::Point3f v1 = (*tempVertexPointer).P() - (*vi).P();
						//vector from gradient field of w
						vcg::Point3f v2 = (*firstVertex).P() - (*vi).P();
						double v1_magnitude = sqrt(v1.X()*v1.X() + v1.Y()*v1.Y() + v1.Z()*v1.Z());
						double v2_magnitude = sqrt(v2.X()*v2.X() + v2.Y()*v2.Y() + v2.Z()*v2.Z());
						v1 = v1 * (1.0 / v1_magnitude);
						v2 = v2 * (1.0 / v2_magnitude);
						double cosTheta = 0.0;
						if (v1_magnitude*v2_magnitude == 0) {
							cosTheta = 0.0;
						}
						else {
							cosTheta = v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
						}
						float tempAngleinDegrees = acos(cosTheta) * 180.0 / PI;
						angle = angle + tempAngleinDegrees;
					}
				}
				if (i > 0) {
					if (previousSet == 'F' && currentSet == 'F') {
						//don't add the angle
					}
					else {
						//calculate and add the angle between previousVertex and tempVertexPointer
						vcg::Point3f v1 = (*tempVertexPointer).P() - (*vi).P();
						//vector from gradient field of w
						vcg::Point3f v2 = (*previousVertex).P() - (*vi).P();
						double v1_magnitude = sqrt(v1.X()*v1.X() + v1.Y()*v1.Y() + v1.Z()*v1.Z());
						double v2_magnitude = sqrt(v2.X()*v2.X() + v2.Y()*v2.Y() + v2.Z()*v2.Z());
						v1 = v1 * (1.0 / v1_magnitude);
						v2 = v2 * (1.0 / v2_magnitude);
						double cosTheta = 0.0;
						if (v1_magnitude*v2_magnitude == 0) {
							cosTheta = 0.0;
						}
						else {
							cosTheta = v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z();
						}
						float tempAngleinDegrees = acos(cosTheta) * 180.0 / PI;
						angle = angle + tempAngleinDegrees;
					}
					//update the pointers
					previousVertex = tempVertexPointer;
					previousSet = setVFB[&(*tempVertexPointer)];
				}
				float tempDivergence = divfval[&(*tempVertexPointer)];
				sumOfDivergence = sumOfDivergence + tempDivergence;
				numberOfNeighbors++;
			}
			//Candidate connecting point
			if (numberOfChanges <= 2) {
				//not element of Fp, may be connecting point
				if (angle > connectingPointAngle) {
					(*vi).C() = vcg::Color4b(vcg::Color4b::Blue);
					CpointOwner[vi] = CPointCounter;
					source[vi] = &(*vi);
					vectCPoints.push_back(&(*vi));
					CPointCounter++;
				}
			}
		}//in set 'F'
		vertexCounter++;
	}//go through vertices

	return vectCPoints;
}

vector<CVertexO *> ToothSegmentationPlugin::FindMidPoints(MeshModel &m, RichParameterSet &par, std::vector<CVertexO*> vectCPoints)
{
	// Find Mid Points
	vector<CMeshO::VertexPointer> midPointVector;

	// ****** Find connecting paths from connecting points
	// setup the seed vector of connecting points
	std::vector<VertDist> seedVec;
	std::vector<CVertexO*>::iterator fi;
	for (fi = vectCPoints.begin(); fi != vectCPoints.end(); ++fi)
	{
		seedVec.push_back(VertDist(*fi, 0.0));
	}
	std::vector<VertDist> frontier;
	CMeshO::VertexPointer curr;
	CMeshO::ScalarType unreached = std::numeric_limits<CMeshO::ScalarType>::max();
	CMeshO::VertexPointer pw;
	TempDataType TD(m.cm.vert, unreached);
	std::vector <VertDist >::iterator ifr;
	for (ifr = seedVec.begin(); ifr != seedVec.end(); ++ifr) {
		TD[(*ifr).v].d = 0.0;
		(*ifr).d = 0.0;
		TD[(*ifr).v].source = (*ifr).v;
		frontier.push_back(VertDist((*ifr).v, 0.0));
	}
	// initialize Heap
	make_heap(frontier.begin(), frontier.end(), pred());
	std::vector<int> doneConnectingPoints;
	CMeshO::ScalarType curr_d, d_curr = 0.0, d_heap;
	CMeshO::VertexPointer curr_s = NULL;
	CMeshO::PerVertexAttributeHandle <CMeshO::VertexPointer> * vertSource = NULL;
	CMeshO::ScalarType max_distance = 0.0;
	CMeshO::ScalarType maxdist = par.getFloat("MaxDist");   // ��������Χ
	std::vector<VertDist >::iterator iv;
	int connectionCounter = 0;
	while (!frontier.empty() && max_distance < maxdist)
	{
		pop_heap(frontier.begin(), frontier.end(), pred());
		curr = (frontier.back()).v;
		int tempCpointOwner = CpointOwner[curr];
		int tempFSetNum = FSetNum[curr];
		curr_s = TD[curr].source;
		if (vertSource != NULL)
			(*vertSource)[curr] = curr_s;
		d_heap = (frontier.back()).d;
		frontier.pop_back();
		assert(TD[curr].d <= d_heap);
		assert(curr_s != NULL);
		if (TD[curr].d < d_heap)// a vertex whose distance has been improved after it was inserted in the queue
			continue;
		assert(TD[curr].d == d_heap);
		d_curr = TD[curr].d;
		if (d_curr > max_distance) {
			max_distance = d_curr;
		}
		//check the vertices around the current point
		(*curr).SetV();
		vcg::tri::UpdateNring<CMeshO> OneRing(&(*curr), &m.cm);
		OneRing.insertAndFlag1Ring(&(*curr));
		(*curr).ClearV();
		for (int i = 0; i < OneRing.allV.size(); i++) {
			pw = OneRing.allV.at(i);
			//just find the shortest distance between points
			curr_d = d_curr + vcg::Distance(pw->cP(), curr->cP());
			//check if we are still searching from this connecting point
			std::vector<int>::iterator it;
			it = std::find(doneConnectingPoints.begin(), doneConnectingPoints.end(), CpointOwner[curr]);
			//we are still searching from this connecting point
			if (it == doneConnectingPoints.end()) {
				//This point has been explored by a different Cpoint before
				//deleted the bit about not connecting to your own set: && FSetNum[pw] != FSetNum[curr]

				if (CpointOwner[pw] != CpointOwner[curr] && CpointOwner[pw] != -1 && setVFB[&(*pw)] != 'F') {
					(*pw).C() = vcg::Color4b(vcg::Color4b::Blue);
					connectionCounter++;
					source1[pw] = curr;
					//store the midpoint so we can make the connecting path later
					midPointVector.push_back(pw);
					//add CpointOwner to vector of doneConnectingPoints
					doneConnectingPoints.push_back(CpointOwner[curr]);
					doneConnectingPoints.push_back(CpointOwner[pw]);
				}
				//deleted && curvature[pw] < maxCurvatureInPath
				else if (TD[(pw)].d > curr_d && curr_d < maxdist && setVFB[&(*pw)] != 'F') {
					//This point has not been explored before, keep looking
					//update source, Fsetnum, CpointOwner
					CpointOwner[pw] = tempCpointOwner;
					source[pw] = curr;
					FSetNum[pw] = tempFSetNum;
					TD[(pw)].d = curr_d;
					TD[pw].source = curr;
					frontier.push_back(VertDist(pw, curr_d));
					push_heap(frontier.begin(), frontier.end(), pred());
				}
			}
			else {
				//not searching from this connecting point anymore
			}
		}
	}// end while
	 //find mid points and color magenta

	return midPointVector;
}

void ToothSegmentationPlugin::AddPathtoConnect(MeshModel &m, std::vector<CVertexO*> midPointVector)
{
	// ****** Make the connecting path
	vector<CMeshO::VertexPointer>::iterator iter;
	for (iter = midPointVector.begin(); iter != midPointVector.end(); ++iter) {
		//track back source to beginning
		CMeshO::VertexPointer tempVertexPointer = *iter;
		CMeshO::VertexPointer tempVertexSource = source[tempVertexPointer];
		while (tempVertexPointer != tempVertexSource) {
			(*tempVertexPointer).C() = vcg::Color4b(vcg::Color4b::Red);
			setVFB[tempVertexPointer] = 'F';
			tempVertexPointer = tempVertexSource;
			tempVertexSource = source[tempVertexPointer];
		}
		//track back source1 to beginning
		tempVertexPointer = *iter;
		tempVertexSource = source1[tempVertexPointer];
		while (tempVertexPointer != tempVertexSource) {
			(*tempVertexPointer).C() = vcg::Color4b(vcg::Color4b::Cyan);
			setVFB[tempVertexPointer] = 'F';
			tempVertexPointer = tempVertexSource;
			tempVertexSource = source[tempVertexPointer];
		}
	}//make connecting path

	 // ����� 'F' �ı��
	CMeshO::VertexIterator vi;
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//if it is not a feature
		if (setVFB[vi] != 'F') {
			setVFB[vi] = 'V';
		}
	}
}

void ToothSegmentationPlugin::RemoveLines(MeshModel &m)
{
	CMeshO::VertexIterator vi;
	//Check if color same on both sides of F line
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (setVFB[vi] == 'F') {
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));
			(*vi).ClearV();
			vcg::Color4b color1, color2, color3;
			vcg::Color4b whiteColor = vcg::Color4b(vcg::Color4b::LightGray);
			vcg::Color4b blueColor = vcg::Color4b(vcg::Color4b::Red);
			vcg::Color4b tempColor;
			int colorCounter = 0;
			CVertexO* tempVertexPointer;
			for (int i = 0; i < OneRing.allV.size(); i++) {
				if (i == 0) {
					tempVertexPointer = OneRing.allV.at(i);
					color1 = (*tempVertexPointer).C();
					colorCounter++;
				}
				if (i > 0 && i < OneRing.allV.size()) {
					tempVertexPointer = OneRing.allV.at(i);
					tempColor = (*tempVertexPointer).C();
					if (colorCounter == 1) {
						if (tempColor != color1) {
							color2 = tempColor;
							colorCounter++;
						}
					}
					if (colorCounter == 2) {
						if (tempColor != color1 && tempColor != color2) {
							color3 = tempColor;
							colorCounter++;
						}
					}
				}
			}
			//delete it from feature
			if (colorCounter == 2 || colorCounter == 1) {
				if (colorCounter == 2) {
					if (color1 != whiteColor && color2 != whiteColor) {
						//delete from feature if surrounded by a color besides white
						setVFB[vi] = 'V';
						if (color1 != blueColor) {
							(*vi).C() = color1;
						}
						if (color2 != blueColor) {
							(*vi).C() = color2;
						}
					}
				}
				if (colorCounter == 1) {
					if (color1 != whiteColor) {
						//delete from feature if surrounded by a color besides white
						setVFB[vi] = 'V';
						if (color1 != blueColor) {
							(*vi).C() = color1;
						}
					}
				}
			}
		}
	}//check if color on both sides of line same
}

MESHLAB_PLUGIN_NAME_EXPORTER(ToothSegmentationPlugin)