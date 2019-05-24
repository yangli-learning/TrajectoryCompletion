#include <iostream>
#include <QSettings>

#include "singleton.h"
#include "parameters.h"

Parameters* Parameters::getInstance(void) {
  return Singleton<Parameters>::instance();
}

Parameters::Parameters(void) :
    q_settings_(NULL) {
  return;
}

Parameters::~Parameters(void) {
  delete q_settings_;

  return;
}

std::string Parameters::getFilename(void) const {
  std::string filename;
  if (q_settings_ != NULL)
    filename = q_settings_->fileName().toStdString();

  return filename;
}

bool Parameters::readParametersImpl(void) {
  bool all_success = true;
#define X(type, name) \
  if (!readParameter(std::string(#name), name))\
  {\
  all_success = false;\
  std::cout << #name << " is uninitialized" << std::endl;\
  name = type();\
  }
  X_GLOBAL_APP_STATE_FIELDS
#undef X

  return all_success;
}

bool Parameters::setParameterFile(const std::string & filename) {
  delete q_settings_;
  q_settings_ = new QSettings(filename.c_str(), QSettings::IniFormat);

  return true;
}

bool Parameters::readParameters(const std::string& filename) {
  setParameterFile(filename);

  return readParametersImpl();
}

bool Parameters::readParameters(void) {
  q_settings_->sync();

  return readParametersImpl();
}

bool Parameters::readParameterImpl(const std::string& name, QVariant& q_variant) {
  QString key(name.c_str());
  if (!q_settings_->contains(key))
    return false;

  q_variant = q_settings_->value(key);
  return true;
}
