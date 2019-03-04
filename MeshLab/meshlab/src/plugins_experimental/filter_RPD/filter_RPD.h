/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2018                                                \/)\/    *
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

#ifndef FILTER_RPD_H
#define FILTER_RPD_H

#include <QObject>
#include <common/interfaces.h>
#include <vcg/math/disjoint_set.h>

using namespace vcg;
using namespace std;

class FilterRPD : public QObject, public MeshFilterInterface
{
    Q_OBJECT
		MESHLAB_PLUGIN_IID_EXPORTER(MESH_FILTER_INTERFACE_IID)
		Q_INTERFACES(MeshFilterInterface)

		enum RefPlane { REF_CENTER, REF_MIN, REF_ORIG };
public:
	/* naming convention :
	- FP -> Filter Plugin
	- name of the filter separated by _
	*/
	enum {
		FP_RPD_IMPORT,
		FP_RPD_EXPORT,
		FP_RPD_AUTO_RECG,
		FP_RPD_MNUL_ADD,
		FP_RPD_MNUL_DELT,
		FP_RPD_GET_GL
	};

	FilterRPD();
	~FilterRPD() {}

	QString filterName(FilterIDType filter) const;
	QString filterInfo(FilterIDType filter) const;

	FilterClass getClass(QAction *);
	void initParameterSet(QAction *, MeshDocument &md, RichParameterSet & /*parent*/);
	bool applyFilter(QAction *filter, MeshDocument &md, RichParameterSet & /*parent*/, vcg::CallBackPos * cb);
	int postCondition(QAction *filter) const;
	int getPreCondition(QAction *filter) const;	
	int getRequirements(QAction *);
	FILTER_ARITY filterArity(QAction *) const { return SINGLE_MESH; }
public:

private:
	bool IsCurComputed;
	struct Teeth {
		Teeth(CVertexO* x) { node = x; count = 1; mark = -1; }
		CVertexO* node;
		int count;
		int mark;
	};
	void ARInit(MeshModel &, bool );
	void ARSmooth(MeshModel &);
	void ARDif(MeshModel &);
	void ARAbstract(MeshModel &, float );
	void ARAreaDilate(MeshModel &, int );
	void ARAreaClean(MeshModel &, RichParameterSet &);
	int ComputeMark(int a);
	int getTMark(vector<pair<CVertexO*, int>>& Mark, CVertexO * vi);
	void setTMark(vector<pair<CVertexO*, int>>& Mark, CVertexO * vi,int a);
	int countBit(int num);
	void transTMark(CVertexO * vi);
public:
	CMeshO::PerVertexAttributeHandle<float> QCopy;
	///顶点标志向量
	vector< pair<CVertexO *, int> > TMark;
	vector< pair<CVertexO *, int> > glTMark;
	///单牙父节点集
	vector<Teeth> Tparent;

	

public:///
	

};

#endif // FILTER_RPD_H