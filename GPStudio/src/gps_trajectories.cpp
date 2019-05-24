#include <stdio.h>
#include <fstream>
#include <limits>
#include <cstdint>

#include <osg/Image>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Texture2D>

#include <boost/filesystem.hpp>

#include <GeographicLib/GeoCoords.hpp>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>

#include <QDebug>

#include "bbox_2d.h"
#include "parameters.h"
#include "gps_trajectory.pb.h"

#include "gps_trajectories.h"

GPSTrajectories::GPSTrajectories(void) :
    show_traj_points_(true) {
	
}

GPSTrajectories::~GPSTrajectories(void) {
}

static void convertCoords(TrajPoint* traj_point) {
  GeographicLib::GeoCoords geo_coords(traj_point->lat() / 1e5, traj_point->lon() / 1e5);
  traj_point->set_x(geo_coords.Easting());
  traj_point->set_y(geo_coords.Northing());
  return;
}

bool GPSTrajectories::loadPBF(const std::string& filename) {
	uint32_t MAX_TRAJ =  Parameters::getInstance()->max_n_trajectories;//10000;
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  std::ifstream fin(filename.c_str(), std::ios::binary);
  if (!fin.good()) {
    fprintf(stderr, "ERROR! Cannot open protobuf trajectory file %s\n",
        filename.c_str());
    return false;
  }
  google::protobuf::io::ZeroCopyInputStream *raw_input = new google::protobuf::io::IstreamInputStream(&fin);
  google::protobuf::io::CodedInputStream *coded_input = new google::protobuf::io::CodedInputStream(raw_input);
  coded_input->SetTotalBytesLimit(1000000000, 500000000);

  uint32_t n_traj;
  if (!coded_input->ReadLittleEndian32(&n_traj)) {
    return false; // possibly end of file
  }
  n_traj = std::min(n_traj,MAX_TRAJ);
  qDebug()<<"load " << n_traj <<" trajectories out of , max is " << MAX_TRAJ;
  gps_trajs_.clear();
  for (size_t i = 0; i < n_traj; ++i) {
    uint32_t length;

    if (!coded_input->ReadLittleEndian32(&length)) {
      break; // possibly end of file
    }

    auto limit = coded_input->PushLimit(length);
    gps_trajs_.push_back(std::shared_ptr<GpsTraj>(new GpsTraj()));
    if (!gps_trajs_.back()->MergePartialFromCodedStream(coded_input)) {
      fprintf(stderr, "ERROR! Protobuf trajectory file contaminated!\n");
      gps_trajs_.pop_back();
      return false;
    } else {
      std::shared_ptr<GpsTraj> gps_traj = gps_trajs_.back();
      for(int j = 0, j_end = gps_traj->point_size(); j < j_end; ++ j) {
        convertCoords(gps_traj->mutable_point(j));
      }
    }
    coded_input->PopLimit(limit);
  }

  fin.close();
  delete coded_input;
  delete raw_input;
  //google::protobuf::ShutdownProtobufLibrary();

  return true;
}
bool GPSTrajectories::loadTXT(const std::string& filename) {
  FILE* fin = fopen(filename.c_str(), "r");
  if (fin == nullptr) {
    fprintf(stderr, "ERROR! Cannot open txt trajectory file %s\n", filename.c_str());
    return false;
  }
  int n_traj = 0;
  fscanf(fin, "%d", &n_traj);

  for (int i = 0; i < n_traj; ++i) {
    int n_points = 0;
    gps_trajs_.push_back(std::shared_ptr < GpsTraj > (new GpsTraj()));
    std::shared_ptr<GpsTraj> gps_traj = gps_trajs_.back();
    fscanf(fin, "%d", &n_points);
    for (int j = 0; j < n_points; ++j) {
      float x, y;
      fscanf(fin, "%f %f", &x, &y);
      TrajPoint* traj_point = gps_traj->add_point();
      traj_point->set_x(x);
      traj_point->set_y(y);
    }
  }

  return true;
}

void GPSTrajectories::toggleRenderTrajPoints(bool show_traj_points) {
  if (show_traj_points_ == show_traj_points)
    return;

  QWriteLocker lock(&read_write_lock_);
  expired_ = true;

  show_traj_points_ = !show_traj_points_;

  return;
}

bool GPSTrajectories::load(const std::string& filename) {
   QWriteLocker lock(&read_write_lock_);
	//	std::cout <<filename <<"DDDDDD";
  expired_ = true;
 // std::string extension = boost::filesystem::path(filename).extension().string();
  bool success = false;
 // if (extension == ".pbf") {
    success = loadPBF(filename);
  //} else if (extension == ".txt") {
  //  success = loadTXT(filename);
  //}

  if (success) {
    float min_x = std::numeric_limits<float>::max();
    float min_y = std::numeric_limits<float>::max();
    float max_x = std::numeric_limits<float>::min();
    float max_y = std::numeric_limits<float>::min();
 
    for (size_t i = 0, i_end = gps_trajs_.size(); i < i_end; ++i) {
      std::shared_ptr<GpsTraj>& gps_traj = gps_trajs_[i];
      for (size_t j = 0, j_end = gps_traj->point_size(); j < j_end; ++j) {
        float x = gps_traj->point(j).x();
        float y = gps_traj->point(j).y();
        min_x = std::min(min_x, x);
        max_x = std::max(max_x, x);
        min_y = std::min(min_y, y);
        max_y = std::max(max_y, y);
      }
    }
    min_corner_ = osg::Vec2(min_x, min_y);
    max_corner_ = osg::Vec2(max_x, max_y);
  }

  return success;
} 

void GPSTrajectories::updateImpl(void) {
  if (!show_traj_points_)
    return;

  osg::ref_ptr<osg::Vec3Array> vertices = new osg::Vec3Array;
  osg::ref_ptr<osg::Vec3Array> normals = new osg::Vec3Array;
  osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;

  for (size_t i = 0, i_end = gps_trajs_.size(); i < i_end; ++i) {
    std::shared_ptr<GpsTraj>& gps_traj = gps_trajs_[i];
    for (size_t j = 0, j_end = gps_traj->point_size(); j < j_end; ++j) {
      vertices->push_back(osg::Vec3(gps_traj->point(j).x(), gps_traj->point(j).y(), 0.0f));
      normals->push_back(osg::Vec3(0.0f, 0.0f, 1.0f));
      colors->push_back(osg::Vec4(0.8f, 0.8f, 0.8f, 1.0f));//0.8
    }
  }

  if (!vertices->empty()) {
    osg::ref_ptr<osg::Geometry> geometry = new osg::Geometry;
    geometry->setUseDisplayList(true);
    geometry->setVertexArray(vertices);
    geometry->setNormalArray(normals);
    normals->setBinding(osg::Array::BIND_PER_VERTEX);
    geometry->setColorArray(colors);
    colors->setBinding(osg::Array::BIND_PER_VERTEX);
    geometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, vertices->size()));
    osg::Geode* geode = new osg::Geode;
    geode->addDrawable(geometry);
    content_root_->addChild(geode);
  }

  return;
}


