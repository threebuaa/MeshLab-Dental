/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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

#ifndef TOOTH_SEGMENTATIONPLUGIN_H
#define TOOTH_SEGMENTATIONPLUGIN_H

#include <common/interfaces.h>

// �μ� geodesic.h 167����������
/* Auxiliary stcuct for heap of vertices to visit while searching connecting points
*/
struct VertDist {
	VertDist() {}
	VertDist(CMeshO::VertexPointer _v, float _d) :v(_v), d(_d) {}
	CMeshO::VertexPointer v;
	CMeshO::ScalarType d;
};
/* Temporary data to for connecting points search: estimated distance and source
*/
struct TempData {
	TempData() {}
	TempData(const CMeshO::ScalarType & d_) { d = d_; source = NULL; parent = NULL; }
	CMeshO::ScalarType d;
	CMeshO::VertexPointer source;//closest source
	CMeshO::VertexPointer parent;
};

typedef vcg::SimpleTempData< std::vector<CMeshO::VertexType>, TempData > TempDataType;

struct pred : public std::binary_function<VertDist, VertDist, bool> {
	pred() {};
	bool operator()(const VertDist& v0, const VertDist& v1) const
	{
		return (v0.d > v1.d);
	}
};

class ToothSegmentationPlugin : public QObject, public MeshFilterInterface
{
    Q_OBJECT
        MESHLAB_PLUGIN_IID_EXPORTER(MESH_FILTER_INTERFACE_IID)
        Q_INTERFACES(MeshFilterInterface)

        enum RefPlane { REF_CENTER,REF_MIN,REF_ORIG};

public:
    /* naming convention :
    - FP -> Filter Plugin
    - name of the filter separated by _
    */
    enum {
		FP_MESH_SEGMENTATION
    } ;


    ToothSegmentationPlugin();
    ~ToothSegmentationPlugin(){}
	virtual QString pluginName(void) const { return "ToothSengmentationPlugin"; }

    QString filterName(FilterIDType filter) const;
    QString filterInfo(FilterIDType filter) const;

    FilterClass getClass(QAction *);
    void initParameterSet(QAction *, MeshModel & /*m*/, RichParameterSet & /*parent*/);
    bool applyFilter(QAction *filter, MeshDocument &md, RichParameterSet & /*parent*/, vcg::CallBackPos * cb) ;
    int postCondition(QAction *filter) const;
    //int getPreCondition(QAction *filter) const;
    FILTER_ARITY filterArity(QAction *) const {return SINGLE_MESH;}

public:
	bool CurvGroup(MeshModel &m, RichParameterSet &par);    // �����ʾ�����ʴ���ĳһ��ֵ�Ķ���, ���� FeaturePoints ��
	//void DelUnuseSection(MeshModel &m); // ɾ��ģ���϶������õĲ���
	void ShowFeature(MeshModel &m); // ��ʾ��ȡ������

	// ��̬ѧ��������
	void meshDilate(MeshModel &m, int n_ring);  // dilate ����
	void meshErode(MeshModel &m, int n_ring);   // erode ����
	void meshSkeletonize(MeshModel &m); // skeletonize ����
	void meshPrune(MeshModel &m, int pruneLength);   // prune ����

	// ��ȡ
	void FilterNoise(MeshModel &m, float deletePercentOfF);   // ȥ�� Feature Skeleton �ϵ�СƬ��������	
	void showSegments1st(MeshModel &m);
	void CalculateDistanceField(MeshModel &m, std::vector<CVertexO*> vectCPoints);

	std::vector<CVertexO*> FindConnectPoints(MeshModel &m);
	std::vector<CVertexO*> FindMidPoints(MeshModel &m, RichParameterSet &par, std::vector<CVertexO*> vectCPoints);
	void CalculateDistanceField(MeshModel &m);
	void AddPathtoConnect(MeshModel &m, std::vector<CVertexO*> midPointVector);
	void RemoveLines(MeshModel &m); // ȥ�������Ƿ�յ�������

	void SaveAsSTL(MeshModel &m);

public:
	std::vector<CMeshO::VertexPointer> FeaturePoints;
	float CurvThreshold;

public:
	// �Զ��嶥������
	///��������,Vertex, Feature, Border
	CMeshO::PerVertexAttributeHandle<char> setVFB;	
	CMeshO::PerVertexAttributeHandle<int> marked;
	CMeshO::PerVertexAttributeHandle<float> curvature;
	CMeshO::PerVertexAttributeHandle<int> disc;
	CMeshO::PerVertexAttributeHandle<int> center;
	CMeshO::PerVertexAttributeHandle<int> FSetNum;
	CMeshO::PerVertexAttributeHandle<float> disth; // ������볡���������, �洢ÿ�㵽�߽��ߵľ���
	CMeshO::PerVertexAttributeHandle<int> CpointOwner;  // connecting points, ��Ҫ������ĩ�˵�
	CMeshO::PerVertexAttributeHandle<float> divfval;    // "DivergenceFieldValue", �ݶȳ�ֵ
	CMeshO::PerVertexAttributeHandle<CMeshO::VertexPointer> source; // ����ĩ�˵�����·���ϵĵ�
	CMeshO::PerVertexAttributeHandle<CMeshO::VertexPointer> source1; // VertexPointer named "source1"
};

#endif
