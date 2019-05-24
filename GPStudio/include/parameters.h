#pragma once

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <stdexcept>

#include <QVariant>

#define X_GLOBAL_APP_STATE_FIELDS \
  X(bool, dump_query_image) \
  X(float, map_image_resolution) \
  X(int, n_directions) \
  X(bool, heading_linear_interpolation) \
  X(bool, skeleton_use_heading)\
  X(float, skeleton_thresh) \
  X(float, skeleton_mu) \
  X(float, skeleton_h0) \
  X(int, skeleton_sample_cell_size) \
  X(int, skeleton_n_samples_per_cell) \
  X(int, skeleton_max_iter)\
  X(bool, skeleton_normalize_heading)\
  X(float, skeleton_heading_thresh)\
  X(float, skeleton_sigma_thresh)\
  X(float, skeleton_branch_cos_thresh)\
  X(float, skeleton_max_trace_step)\
  X(float, skeleton_noise_scale)\
  X(int, skeleton_max_neighbors)\
  X(float, skeleton_branch_min_angle)\
  X(int, max_n_trajectories)\
  X(bool, heading_blur)\
  X(float, skeleton_cluster_radius)\
  X(bool, skeleton_position_pca)\
  X(bool,skeleton_mvn_kernel)\
  X(int, skeleton_linearity_function)\
  X(float, skeleton_proj_angle_thresh)\
  X(float, skeleton_proj_ratio_thresh)\
  X(float, skeleton_extension_step)\
  X(int, skeleton_extension_max_iter)\
  X(float, skeleton_junction_smoothness)\
  X(float, skeleton_min_edge_len)\
  X(float, skeleton_max_edge_len)\
  X(int, skeleton_min_cluster_size)\
  X(int, cluster_graph_deg)\
  X(float,cluster_graph_ang_thresh)\
  X(float,cluster_graph_neighbor_penalty)\
  X(bool,skeleton_edge_smoothing)\
  X(float,min_flow_ratio)

#ifndef VAR_NAME
#define VAR_NAME(x) #x
#endif

class QSettings;

Q_DECLARE_METATYPE(std::string);

class Parameters {
public:
  Parameters(void);
  virtual ~Parameters(void);

  static Parameters* getInstance(void);

#define X(type, name) type name;
  X_GLOBAL_APP_STATE_FIELDS
#undef X

  std::string getFilename(void) const;
  bool setParameterFile(const std::string & filename);
  bool readParameters(const std::string& filename);
  bool readParameters(void);

protected:
  template<class T>
  bool readParameter(const std::string& name, T& value) {
    QVariant q_variant;
    if (!readParameterImpl(name, q_variant))
      return false;

    if (!q_variant.canConvert<T>()) {
      // try to interpret as std::string...
      q_variant.setValue(q_variant.value<QString>().toStdString());
    }

    if (!q_variant.canConvert<T>())
      return false;

    value = q_variant.value<T>();
    return true;
  }
private:
  bool readParametersImpl(void);
  bool readParameterImpl(const std::string& name, QVariant& q_variant);

  QSettings* q_settings_;
};

#endif // PARAMETERS_H
