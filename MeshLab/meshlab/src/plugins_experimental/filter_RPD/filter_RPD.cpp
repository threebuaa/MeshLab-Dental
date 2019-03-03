/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2018                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
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

#include <filter_RPD.h>
#include <math.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>
#include <vcg/complex/algorithms/update/nring.h>


#include <fstream>
#include <iostream>
#include <iomanip>

using namespace vcg;
using namespace std;


FilterRPD::FilterRPD()
{
	typeList << FP_RPD_IMPORT
		<< FP_RPD_EXPORT
		<< FP_RPD_AUTO_RECG
		<< FP_RPD_MNUL_ADD
		<< FP_RPD_MNUL_DELT
		<< FP_RPD_GET_GL;

	foreach(FilterIDType tt, types())
		actionList << new QAction(filterName(tt), this);
	//test = 0;
}

//FilterRPD::~FilterRPD() {	}

QString FilterRPD::filterInfo(FilterIDType filter) const
{
	switch (filter)
	{
	case FP_RPD_IMPORT:					return QString("Import Model...");
	case FP_RPD_EXPORT:            		return QString("Export Model\nx,y,z,k1,k2,H,G");
	case FP_RPD_AUTO_RECG:				return QString("Intelligent recognition GL_Area according to vertex curvature");
	case FP_RPD_MNUL_ADD:				return QString("Add selected area into GL_Area");
	case FP_RPD_MNUL_DELT:				return QString("Delete selected area from GL_Area");
	case FP_RPD_GET_GL:					return QString("Get Gingival Line by regional growth");

	default: assert(0);
	}
	return QString("error!");
}

QString FilterRPD::filterName(FilterIDType filter) const
{
	switch (filter)
	{
	case FP_RPD_IMPORT:					return tr("Import Model...");
	case FP_RPD_EXPORT:            		return tr("Export Model");
	case FP_RPD_AUTO_RECG:				return tr("Intelligent Recognition GL_Area");
	case FP_RPD_MNUL_ADD:				return tr("Add GL_Area Manually");
	case FP_RPD_MNUL_DELT:				return tr("Delete GL_Area Manually");
	case FP_RPD_GET_GL:					return tr("Get Gingival Line");

	default: assert(0);
	}
	return QString("error!");
}

MeshFilterInterface::FilterClass FilterRPD::getClass(QAction* action)
{
	switch (ID(action))
	{
	case FP_RPD_IMPORT:
	case FP_RPD_EXPORT:
		return MeshFilterInterface::RPDDebug;
	case FP_RPD_AUTO_RECG:
	case FP_RPD_MNUL_ADD:
	case FP_RPD_MNUL_DELT:
	case FP_RPD_GET_GL:
		return MeshFilterInterface::RPDGL;
	}	
}

int FilterRPD::postCondition(QAction * filter) const
{
	switch (ID(filter))
	{
	case FP_RPD_AUTO_RECG:
	case FP_RPD_MNUL_ADD:
	case FP_RPD_MNUL_DELT:
	case FP_RPD_GET_GL:
		return MeshModel::MM_VERTCOLOR;
	default: return 0;
	}
}

int FilterRPD::getPreCondition(QAction * filter) const
{
	return 0;
}

int FilterRPD::getRequirements(QAction * filter)
{
	switch (ID(filter))
	{
	case FP_RPD_AUTO_RECG:
	case FP_RPD_MNUL_ADD:
	case FP_RPD_MNUL_DELT:
	case FP_RPD_GET_GL:
		return MeshModel::MM_VERTCOLOR;
	default: return 0;
	}
}


void FilterRPD::initParameterSet(QAction *action, MeshDocument &md, RichParameterSet & parlst)
{
	int id = ID(action);
	if (id == FP_RPD_AUTO_RECG)
	{
		parlst.addParam(new RichBool("Selected", false, "Selected only", "If checked, only recognize vertex in selected areas."));

		QStringList lst;
		lst << "K2"<< "Maen";
		parlst.addParam(new RichEnum("CurvatureType", 0, lst, "Curvature type",	"The type of the curvature to plot."));
		
		parlst.addParam(new RichInt("MinAreaP", 200, "Threshold", "Threshold of the Minimum Area Points"));
		return;
	}
/*	if (id == FP_MORPHOLOGY_OP)
	{
		parlst.addParam(new RichInt("DilateCount", 2, "Dilate count", "The number of model dilation operations."));
		parlst.addParam(new RichInt("ErodeCount", 1, "Erode count", "The number of model erosion operations."));
	}*/
}

bool FilterRPD::applyFilter(QAction *filter, MeshDocument &md, RichParameterSet & par, vcg::CallBackPos * cb)
{
	if (md.mm() == NULL)	return false;
	
	MeshModel &m = *(md.mm());
	
	CMeshO::VertexIterator vi;

	switch (ID(filter))
	{
	case FP_RPD_MNUL_ADD:
	{
		tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
		tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(m.cm);	
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++) {
			if (vi->IsS() == true)
			{
				vi->SetGL();
				vi->C() = vcg::Color4b::Red;
			}
		}
		tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
		tri::UpdateSelection<CMeshO>::FaceClear(m.cm);
	}break;
	case FP_RPD_MNUL_DELT:
	{
		tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
		tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(m.cm);		
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++) {
			if (vi->IsS() == true)
			{
				vi->ClearGL();
				vi->C() = vcg::Color4b::White;
			}
		}
		tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
		tri::UpdateSelection<CMeshO>::FaceClear(m.cm);
	}break;
	case FP_RPD_AUTO_RECG:
	{	
		for (int i = 0; i < md.size(); i++)
		{
			MeshModel &m = *md.meshList.at(i);
			if (!m.isVisible())
				continue;
			//filter准备
			m.updateDataMask(MeshModel::MM_VERTCURVDIR);
			m.updateDataMask(MeshModel::MM_VERTFACETOPO | MeshModel::MM_FACEFACETOPO);
			m.updateDataMask(MeshModel::MM_VERTQUALITY | MeshModel::MM_VERTCOLOR);
			QCopy = vcg::tri::Allocator<CMeshO>::AddPerVertexAttribute<float>(m.cm, string("QCopy"));
			//流形检查
			if (vcg::tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m.cm) > 0 || vcg::tri::Clean<CMeshO>::CountNonManifoldVertexFF(m.cm) > 0)
			{
				this->errorMessage = "Mesh has some not 2-manifold faces, cannot compute principal curvature directions";
				return false; // can't continue, mesh can't be processed
			}
			//同步交互选择结果
			tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
			tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(m.cm);
			//曲率计算
			vcg::tri::UpdateCurvatureFitting<CMeshO>::computeCurvature(m.cm);
			//确定曲率类型
			switch (par.getEnum("CurvatureType"))
			{
			case 0:		vcg::tri::UpdateQuality<CMeshO>::VertexFromK2CurvatureDir(m.cm);	break;
			case 1:		vcg::tri::UpdateQuality<CMeshO>::VertexFromMeanCurvatureDir(m.cm);	break;
			default: assert(0);
			}
			//核心
			if (m.hasDataMask(MeshModel::MM_VERTQUALITY))
			{
				ARInit(m, par.getBool("Selected"));
				ARSmooth(m);
				ARDif(m);
				if (par.getEnum("CurvatureType"))
					ARAbstract(m, -0.50f);
				else
					ARAbstract(m, -1.00f);
				ARAreaDilate(m, 3);
				ARAreaClean(m, par);
			}
			else
			{
				this->errorMessage = "The Mesh has no curvature in vertex quality";
				return false;
			}

			vcg::tri::Allocator<CMeshO>::DeletePerVertexAttribute<float>(m.cm, QCopy);
		}
	}break;
	case FP_RPD_GET_GL:
	{	//filter准备
		m.updateDataMask(MeshModel::MM_VERTMARK);
		m.updateDataMask(MeshModel::MM_VERTFACETOPO | MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		
		CMeshO::VertexIterator vi;		
		//同步GLarea修正结果
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++) {
			Color4<unsigned char> cc = vi->C();
			if ((unsigned int)cc[0] == 255 && (unsigned int)cc[1] == 0 && (unsigned int)cc[2] == 0)
				vi->SetGL();
			else
				vi->ClearGL();
		}	
		//清理
		{
			//waiting...


		}
		//建立并查集
		DisjointSet<CVertexO> *TeethDS = new vcg::DisjointSet<CVertexO>();
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (!vi->IsGL()) {
				TeethDS->MakeSet(&(*vi));   // 初始时, 每个顶点单独作为一个 set
			}
		}//init		
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			//if it is a feature, then find 1-ring and use Union to merge sets
			if (!vi->IsGL()) {
				(*vi).SetV();
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
				OneRing.insertAndFlag1Ring(&(*vi));
				(*vi).ClearV();
				for (int i = 0; i < OneRing.allV.size(); i++) {
					CVertexO* tempVP = OneRing.allV.at(i);
					//if the neighbor of the current vertex is a feature then merge sets
					if (!tempVP->IsGL()) {
						//if they are not already in the same set, merge them
						if (TeethDS->FindSet(&(*vi)) != TeethDS->FindSet(tempVP)) {
							TeethDS->Union(&(*vi), tempVP);
						}
					}
				}
			}
		}//union		
		//count
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (!vi->IsGL()) {
				CVertexO* cRoot = TeethDS->FindSet(&(*vi));
				if (Tparent.empty())
					Tparent.push_back(cRoot);
				else
				{
					int mark = 0;
					for (int i = 0; i < Tparent.size(); i++) {
						if (cRoot == Tparent[i].node)
						{
							Tparent[i].count++;
							mark = 1;
							break;
						}
					}
					if (mark == 0) {
						Tparent.push_back(cRoot);
					}
				}
			}
		}
		//计算mark
		int tcount = 0;
		int tMax = Tparent[0].count;
		Tparent[0].mark = 1;
		for (int i = 1; i < Tparent.size(); i++) {
			if (Tparent[i].count > tMax)
			{
				Tparent[tcount].mark = ComputeMark(i);
				tcount = i;
				tMax = Tparent[i].count;
				Tparent[i].mark = 1;
			}
			else
			{
				Tparent[i].mark = ComputeMark(i);
			}
		}	
		//标志
		for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
			if (!vi->IsGL()) {
				int cmark;
				CVertexO* cRoot = TeethDS->FindSet(&(*vi));
				for (int i = 0; i < Tparent.size(); i++) {
					if (Tparent[i].node == cRoot)
					{
						cmark = Tparent[i].mark;
						break;
					}
				}
				if (cmark > 0)
				{
					int a = 0;
				}
				TMark.push_back(pair<CVertexO*, int>(&(*vi), cmark));
			}
			else
			{
				TMark.push_back(pair<CVertexO*, int>(&(*vi), 0));
			}
		}
		//线化
		int change = 0;
		vector<CVertexO *> delglp;
		do {
			change = 0;
			delglp.clear();
			for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
				if (vi->IsGL()) {
					int tpmark = 0;
					vi->SetV();
					vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
					OneRing.insertAndFlag1Ring(&(*vi));
					vi->ClearV();
					for (int i = 0; i < OneRing.allV.size(); i++) {
						if (!OneRing.allV.at(i)->IsGL())
							tpmark |= getTMark(OneRing.allV.at(i));
					}
					setTMark(&(*vi), tpmark);
					if (countBit(tpmark) == 1) {
						delglp.push_back(&(*vi));
					}
				}
			}
			for (int i = 0; i < delglp.size(); i++) {
				int tpmark = 0;
				delglp[i]->SetV();
				vcg::tri::UpdateNring<CMeshO> OneRing(&(*delglp[i]), &m.cm);
				OneRing.insertAndFlag1Ring(&(*delglp[i]));
				delglp[i]->ClearV();
				for (int i = 0; i < OneRing.allV.size(); i++) {
					tpmark |= getTMark(OneRing.allV.at(i));
				}
				if (tpmark == getTMark(delglp[i])) {
					delglp[i]->ClearGL();
					delglp[i]->C() = vcg::Color4b::White;
					change = 1;
				}				
			}
		} while (change != 0);

	}break;
	case FP_RPD_EXPORT:
	{
		for (int i = 0; i < md.size(); i++)
		{
			md.setCurrentMesh(i);
			MeshModel &m = *(md.mm());
			QString suffix = m.suffixName();
			QString txt = m.label().remove("." + suffix) + ".txt";
			string str = txt.toStdString();
			const char* filename = str.c_str();
			ofstream OutFile(filename);
			OutFile << "point.x	|	point.y	|	point.z	|	k1	|	k2	|	H	|	G	|	ColorR	|	ColorG	|	ColorB\n";
			//filter准备
			m.updateDataMask(MeshModel::MM_VERTCOLOR);
			m.updateDataMask(MeshModel::MM_VERTCURVDIR);
			m.updateDataMask(MeshModel::MM_VERTFACETOPO | MeshModel::MM_FACEFACETOPO);
			//曲率计算
			vcg::tri::UpdateCurvatureFitting<CMeshO>::computeCurvature(m.cm);
			for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
			{
				Color4<unsigned char> cc = vi->C();
				OutFile << setiosflags(ios::fixed | ios::showpoint) << vi->P().X() << "	|	" << vi->P().Y() << "	|	" << vi->P().Z() << "	|	" << vi->K1() << "	|	" << vi->K2() << "	|	" << (vi->K1() + vi->K2()) / 2 << "	|	" << vi->K1()*vi->K2() << "	|	" << (unsigned int)cc[0] << "	|	" << (unsigned int)cc[1] << "	|	" << (unsigned int)cc[2] << endl;
			}
			OutFile.close();
		}
	}break;
	default:  assert(0);
	}
	return true;
}

//粗加工
void FilterRPD::ARInit(MeshModel & m, bool onlyselected)
{
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		vi->C() = vcg::Color4b::White;

		if ((vi->K2() < 0.00f) && !onlyselected)
		{
			vi->SetGL();
			vi->C() = vcg::Color4b::Red;
		}
		if ((vi->K2() < 0.00f) && onlyselected && vi->IsS())
		{
			vi->SetGL();
			vi->C() = vcg::Color4b::Red;
		}
	}
}

//平滑顶点曲率_中值滤波
void FilterRPD::ARSmooth(MeshModel & m)
{
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		if (vi->IsGL())
		{
			vcg::tri::UpdateNring<CMeshO> neighbor(&(*vi), &m.cm);
			neighbor.insertAndFlag1Ring(&(*vi));
			neighbor.expand();

			///中值平滑滤波
			vector<float> tempara;
			for (int i = 0; i < neighbor.allV.size(); i++)
				tempara.push_back(neighbor.allV.at(i)->Q());
			///泡排
			float temp;
			int i, flag = 1;
			i = neighbor.allV.size() - 2;
			while (i > 0 && flag == 1)
			{
				flag = 0;
				for (int j = 0; j <= i; j++)
					if (tempara[j] > tempara[j + 1])
					{
						temp = tempara[j];
						tempara[j] = tempara[j + 1];
						tempara[j + 1] = temp;
						flag = 1;
					}
				i--;
			}///泡排End
			if (neighbor.allV.size() % 2 == 0)
				QCopy[vi] = (tempara[neighbor.allV.size() / 2] + tempara[neighbor.allV.size() / 2 - 1]) / 2;
			else
				QCopy[vi] = tempara[(neighbor.allV.size() - 1) / 2];
		}
	}
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		if (vi->IsGL())
			vi->Q() = QCopy[vi];
		else
			vi->Q() = 0;
	}
}

//增强顶点差异性_拉普拉斯
void FilterRPD::ARDif(MeshModel & m)
{
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		if (vi->IsGL())
		{
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> neighbor(&(*vi), &m.cm);
			neighbor.insertAndFlag1Ring(&(*vi));
			neighbor.expand();
			(*vi).ClearV();

			QCopy[vi] = 0.00f;
			for (int i = 0; i < neighbor.allV.size(); i++)
				QCopy[vi] += neighbor.allV.at(i)->Q();
			QCopy[vi] = QCopy[vi] / neighbor.allV.size() - vi->Q();

			///系数加成
			vector<float> tempara;
			for (int i = 0; i < neighbor.allV.size(); i++)
				tempara.push_back(neighbor.allV.at(i)->Q());
			///up泡排
			float temp;
			int i, flag = 1;
			i = neighbor.allV.size() - 2;
			while (i > 0 && flag == 1)
			{
				flag = 0;
				for (int j = 0; j <= i; j++)
					if (tempara[j] > tempara[j + 1])
					{
						temp = tempara[j];
						tempara[j] = tempara[j + 1];
						tempara[j + 1] = temp;
						flag = 1;
					}
				i--;
			}///泡排End
			if (vi->Q() < tempara[0])
				QCopy[vi] *= 1.5;
			else if (vi->Q() < tempara[1])
				QCopy[vi] *= 1.3;
			else if (vi->Q() < tempara[2])
				QCopy[vi] *= 1.1;
		}
	}
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi)
	{
		if (vi->IsGL())
			vi->Q() -= QCopy[vi];
		else
			vi->Q() = 0.00f;
	}
}

//提取特征区域
void FilterRPD::ARAbstract(MeshModel & m, float f)
{
	for (CMeshO::VertexIterator vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++) {
		if (vi->Q() < f)
		{
			vi->SetGL();
			vi->C() = vcg::Color4b::Red;
		}
		else
		{
			vi->ClearGL();
			vi->C() = vcg::Color4b::White;
		}
	}
}

//Area++
void FilterRPD::ARAreaDilate(MeshModel & m, int count)
{
	CMeshO::VertexIterator vi;
	tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++) {
		if (vi->IsGL())	vi->SetS();
	}
	for (int k = 0; k < count; k++)
	{		
		tri::UpdateSelection<CMeshO>::FaceFromVertexLoose(m.cm);
		tri::UpdateSelection<CMeshO>::VertexFromFaceLoose(m.cm);		
	}
	for (int k = 0; k < count-1; k++)
	{
		tri::UpdateSelection<CMeshO>::VertexFromFaceStrict(m.cm);
		tri::UpdateSelection<CMeshO>::FaceFromVertexStrict(m.cm);
	}
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); vi++) {
		if (vi->IsS())
		{
			vi->SetGL();
			vi->C() = vcg::Color4b::Red;
		}
	}
	tri::UpdateSelection<CMeshO>::VertexClear(m.cm);
	tri::UpdateSelection<CMeshO>::FaceClear(m.cm);
}

//清理区域
void FilterRPD::ARAreaClean(MeshModel & m , RichParameterSet & par)
{
	CMeshO::VertexIterator vi;
	vcg::DisjointSet<CVertexO>* glDJS = new vcg::DisjointSet<CVertexO>();
	vector< pair<CVertexO*, int> > glarea;
	vector<CVertexO*> delt;
	//init
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (vi->IsGL()) {
			glDJS->MakeSet(&(*vi));   // 初始时, 每个顶点单独作为一个 set
		}
	}
	//union
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		//if it is a feature, then find 1-ring and use Union to merge sets
		if (vi->IsGL()) {
			(*vi).SetV();
			vcg::tri::UpdateNring<CMeshO> OneRing(&(*vi), &m.cm);
			OneRing.insertAndFlag1Ring(&(*vi));
			(*vi).ClearV();
			for (int i = 0; i < OneRing.allV.size(); i++) {
				CVertexO* tempVP = OneRing.allV.at(i);
				//if the neighbor of the current vertex is a feature then merge sets
				if (tempVP->IsGL()) {
					//if they are not already in the same set, merge them
					if (glDJS->FindSet(&(*vi)) != glDJS->FindSet(tempVP)) {
						glDJS->Union(&(*vi), tempVP);
					}
				}
			}
		}
	}
	//count
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (vi->IsGL()) {
			CVertexO* cRoot = glDJS->FindSet(&(*vi));
			if (glarea.empty())
				glarea.push_back(pair<CVertexO*, int>(cRoot, 1));
			else
			{
				int mark = 0;
				for (int i = 0; i < glarea.size(); i++) {
					if (cRoot == glarea[i].first)
					{
						glarea[i].second++;
						mark = 1;
						break;
					}
				}
				if (mark == 0) {
					glarea.push_back(pair<CVertexO*, int>(cRoot, 1));
				}
			}
		}
	}
	//delt
	int cMinAreaP = par.getInt("MinAreaP");
	for (int i = 0; i < glarea.size(); i++) {
		if (glarea[i].second < cMinAreaP)	delt.push_back(glarea[i].first);
	}
	for (vi = m.cm.vert.begin(); vi != m.cm.vert.end(); ++vi) {
		if (vi->IsGL()) {
			CVertexO* cRoot = glDJS->FindSet(&(*vi));
			for (int i = 0; i < delt.size(); i++) {
				if (cRoot == delt[i])
				{
					vi->ClearGL();
					vi->C() = vcg::Color4b::White;
				}
			}
		}
	}
}

int FilterRPD::ComputeMark(int a) {
	return{ 1 << a };
}

int FilterRPD::getTMark(CVertexO *vi) {
	for (int i = 0; i < TMark.size(); i++) {
		if (vi == TMark[i].first)
			return TMark[i].second;
	}
}

void FilterRPD::setTMark(CVertexO *vi,int a) {
	for (int i = 0; i < TMark.size(); i++) {
		if (vi == TMark[i].first)
			TMark[i].second = a;
	}
}

int FilterRPD::countBit(int num) {
	int count = 0;
	while (num != 0) {
		num = num & (num - 1);
		count++;
	}
	return count;
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterRPD)