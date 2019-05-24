// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: gps_trajectory.proto

#define INTERNAL_SUPPRESS_PROTOBUF_FIELD_DEPRECATION
#include "gps_trajectory.pb.h"

#include <algorithm>

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/stubs/once.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/wire_format_lite_inl.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)

namespace {

const ::google::protobuf::Descriptor* GpsTraj_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  GpsTraj_reflection_ = NULL;
const ::google::protobuf::Descriptor* TrajPoint_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  TrajPoint_reflection_ = NULL;

}  // namespace


void protobuf_AssignDesc_gps_5ftrajectory_2eproto() {
  protobuf_AddDesc_gps_5ftrajectory_2eproto();
  const ::google::protobuf::FileDescriptor* file =
    ::google::protobuf::DescriptorPool::generated_pool()->FindFileByName(
      "gps_trajectory.proto");
  GOOGLE_CHECK(file != NULL);
  GpsTraj_descriptor_ = file->message_type(0);
  static const int GpsTraj_offsets_[1] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(GpsTraj, point_),
  };
  GpsTraj_reflection_ =
    new ::google::protobuf::internal::GeneratedMessageReflection(
      GpsTraj_descriptor_,
      GpsTraj::default_instance_,
      GpsTraj_offsets_,
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(GpsTraj, _has_bits_[0]),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(GpsTraj, _unknown_fields_),
      -1,
      ::google::protobuf::DescriptorPool::generated_pool(),
      ::google::protobuf::MessageFactory::generated_factory(),
      sizeof(GpsTraj));
  TrajPoint_descriptor_ = file->message_type(1);
  static const int TrajPoint_offsets_[8] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, car_id_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, timestamp_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, lon_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, lat_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, head_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, speed_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, x_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, y_),
  };
  TrajPoint_reflection_ =
    new ::google::protobuf::internal::GeneratedMessageReflection(
      TrajPoint_descriptor_,
      TrajPoint::default_instance_,
      TrajPoint_offsets_,
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, _has_bits_[0]),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(TrajPoint, _unknown_fields_),
      -1,
      ::google::protobuf::DescriptorPool::generated_pool(),
      ::google::protobuf::MessageFactory::generated_factory(),
      sizeof(TrajPoint));
}

namespace {

GOOGLE_PROTOBUF_DECLARE_ONCE(protobuf_AssignDescriptors_once_);
inline void protobuf_AssignDescriptorsOnce() {
  ::google::protobuf::GoogleOnceInit(&protobuf_AssignDescriptors_once_,
                 &protobuf_AssignDesc_gps_5ftrajectory_2eproto);
}

void protobuf_RegisterTypes(const ::std::string&) {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
    GpsTraj_descriptor_, &GpsTraj::default_instance());
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
    TrajPoint_descriptor_, &TrajPoint::default_instance());
}

}  // namespace

void protobuf_ShutdownFile_gps_5ftrajectory_2eproto() {
  delete GpsTraj::default_instance_;
  delete GpsTraj_reflection_;
  delete TrajPoint::default_instance_;
  delete TrajPoint_reflection_;
}

void protobuf_AddDesc_gps_5ftrajectory_2eproto() {
  static bool already_here = false;
  if (already_here) return;
  already_here = true;
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
    "\n\024gps_trajectory.proto\"$\n\007GpsTraj\022\031\n\005poi"
    "nt\030\001 \003(\0132\n.TrajPoint\"{\n\tTrajPoint\022\016\n\006car"
    "_id\030\001 \002(\005\022\021\n\ttimestamp\030\002 \002(\r\022\013\n\003lon\030\003 \002("
    "\005\022\013\n\003lat\030\004 \002(\005\022\014\n\004head\030\005 \002(\005\022\r\n\005speed\030\006 "
    "\002(\005\022\t\n\001x\030\007 \002(\002\022\t\n\001y\030\010 \002(\002", 185);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "gps_trajectory.proto", &protobuf_RegisterTypes);
  GpsTraj::default_instance_ = new GpsTraj();
  TrajPoint::default_instance_ = new TrajPoint();
  GpsTraj::default_instance_->InitAsDefaultInstance();
  TrajPoint::default_instance_->InitAsDefaultInstance();
  ::google::protobuf::internal::OnShutdown(&protobuf_ShutdownFile_gps_5ftrajectory_2eproto);
}

// Force AddDescriptors() to be called at static initialization time.
struct StaticDescriptorInitializer_gps_5ftrajectory_2eproto {
  StaticDescriptorInitializer_gps_5ftrajectory_2eproto() {
    protobuf_AddDesc_gps_5ftrajectory_2eproto();
  }
} static_descriptor_initializer_gps_5ftrajectory_2eproto_;

// ===================================================================

#ifndef _MSC_VER
const int GpsTraj::kPointFieldNumber;
#endif  // !_MSC_VER

GpsTraj::GpsTraj()
  : ::google::protobuf::Message() {
  SharedCtor();
  // @@protoc_insertion_point(constructor:GpsTraj)
}

void GpsTraj::InitAsDefaultInstance() {
}

GpsTraj::GpsTraj(const GpsTraj& from)
  : ::google::protobuf::Message() {
  SharedCtor();
  MergeFrom(from);
  // @@protoc_insertion_point(copy_constructor:GpsTraj)
}

void GpsTraj::SharedCtor() {
  _cached_size_ = 0;
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
}

GpsTraj::~GpsTraj() {
  // @@protoc_insertion_point(destructor:GpsTraj)
  SharedDtor();
}

void GpsTraj::SharedDtor() {
  if (this != default_instance_) {
  }
}

void GpsTraj::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* GpsTraj::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return GpsTraj_descriptor_;
}

const GpsTraj& GpsTraj::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_gps_5ftrajectory_2eproto();
  return *default_instance_;
}

GpsTraj* GpsTraj::default_instance_ = NULL;

GpsTraj* GpsTraj::New() const {
  return new GpsTraj;
}

void GpsTraj::Clear() {
  point_.Clear();
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
  mutable_unknown_fields()->Clear();
}

bool GpsTraj::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:GpsTraj)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoff(127);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // repeated .TrajPoint point = 1;
      case 1: {
        if (tag == 10) {
         parse_point:
          DO_(::google::protobuf::internal::WireFormatLite::ReadMessageNoVirtual(
                input, add_point()));
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(10)) goto parse_point;
        if (input->ExpectAtEnd()) goto success;
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0 ||
            ::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:GpsTraj)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:GpsTraj)
  return false;
#undef DO_
}

void GpsTraj::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:GpsTraj)
  // repeated .TrajPoint point = 1;
  for (int i = 0; i < this->point_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteMessageMaybeToArray(
      1, this->point(i), output);
  }

  if (!unknown_fields().empty()) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        unknown_fields(), output);
  }
  // @@protoc_insertion_point(serialize_end:GpsTraj)
}

::google::protobuf::uint8* GpsTraj::SerializeWithCachedSizesToArray(
    ::google::protobuf::uint8* target) const {
  // @@protoc_insertion_point(serialize_to_array_start:GpsTraj)
  // repeated .TrajPoint point = 1;
  for (int i = 0; i < this->point_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteMessageNoVirtualToArray(
        1, this->point(i), target);
  }

  if (!unknown_fields().empty()) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        unknown_fields(), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:GpsTraj)
  return target;
}

int GpsTraj::ByteSize() const {
  int total_size = 0;

  // repeated .TrajPoint point = 1;
  total_size += 1 * this->point_size();
  for (int i = 0; i < this->point_size(); i++) {
    total_size +=
      ::google::protobuf::internal::WireFormatLite::MessageSizeNoVirtual(
        this->point(i));
  }

  if (!unknown_fields().empty()) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        unknown_fields());
  }
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void GpsTraj::MergeFrom(const ::google::protobuf::Message& from) {
  GOOGLE_CHECK_NE(&from, this);
  const GpsTraj* source =
    ::google::protobuf::internal::dynamic_cast_if_available<const GpsTraj*>(
      &from);
  if (source == NULL) {
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
    MergeFrom(*source);
  }
}

void GpsTraj::MergeFrom(const GpsTraj& from) {
  GOOGLE_CHECK_NE(&from, this);
  point_.MergeFrom(from.point_);
  mutable_unknown_fields()->MergeFrom(from.unknown_fields());
}

void GpsTraj::CopyFrom(const ::google::protobuf::Message& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void GpsTraj::CopyFrom(const GpsTraj& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool GpsTraj::IsInitialized() const {

  if (!::google::protobuf::internal::AllAreInitialized(this->point())) return false;
  return true;
}

void GpsTraj::Swap(GpsTraj* other) {
  if (other != this) {
    point_.Swap(&other->point_);
    std::swap(_has_bits_[0], other->_has_bits_[0]);
    _unknown_fields_.Swap(&other->_unknown_fields_);
    std::swap(_cached_size_, other->_cached_size_);
  }
}

::google::protobuf::Metadata GpsTraj::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = GpsTraj_descriptor_;
  metadata.reflection = GpsTraj_reflection_;
  return metadata;
}


// ===================================================================

#ifndef _MSC_VER
const int TrajPoint::kCarIdFieldNumber;
const int TrajPoint::kTimestampFieldNumber;
const int TrajPoint::kLonFieldNumber;
const int TrajPoint::kLatFieldNumber;
const int TrajPoint::kHeadFieldNumber;
const int TrajPoint::kSpeedFieldNumber;
const int TrajPoint::kXFieldNumber;
const int TrajPoint::kYFieldNumber;
#endif  // !_MSC_VER

TrajPoint::TrajPoint()
  : ::google::protobuf::Message() {
  SharedCtor();
  // @@protoc_insertion_point(constructor:TrajPoint)
}

void TrajPoint::InitAsDefaultInstance() {
}

TrajPoint::TrajPoint(const TrajPoint& from)
  : ::google::protobuf::Message() {
  SharedCtor();
  MergeFrom(from);
  // @@protoc_insertion_point(copy_constructor:TrajPoint)
}

void TrajPoint::SharedCtor() {
  _cached_size_ = 0;
  car_id_ = 0;
  timestamp_ = 0u;
  lon_ = 0;
  lat_ = 0;
  head_ = 0;
  speed_ = 0;
  x_ = 0;
  y_ = 0;
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
}

TrajPoint::~TrajPoint() {
  // @@protoc_insertion_point(destructor:TrajPoint)
  SharedDtor();
}

void TrajPoint::SharedDtor() {
  if (this != default_instance_) {
  }
}

void TrajPoint::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* TrajPoint::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return TrajPoint_descriptor_;
}

const TrajPoint& TrajPoint::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_gps_5ftrajectory_2eproto();
  return *default_instance_;
}

TrajPoint* TrajPoint::default_instance_ = NULL;

TrajPoint* TrajPoint::New() const {
  return new TrajPoint;
}

void TrajPoint::Clear() {
#define OFFSET_OF_FIELD_(f) (reinterpret_cast<char*>(      \
  &reinterpret_cast<TrajPoint*>(16)->f) - \
   reinterpret_cast<char*>(16))

#define ZR_(first, last) do {                              \
    size_t f = OFFSET_OF_FIELD_(first);                    \
    size_t n = OFFSET_OF_FIELD_(last) - f + sizeof(last);  \
    ::memset(&first, 0, n);                                \
  } while (0)

  if (_has_bits_[0 / 32] & 255) {
    ZR_(car_id_, y_);
  }

#undef OFFSET_OF_FIELD_
#undef ZR_

  ::memset(_has_bits_, 0, sizeof(_has_bits_));
  mutable_unknown_fields()->Clear();
}

bool TrajPoint::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!(EXPRESSION)) goto failure
  ::google::protobuf::uint32 tag;
  // @@protoc_insertion_point(parse_start:TrajPoint)
  for (;;) {
    ::std::pair< ::google::protobuf::uint32, bool> p = input->ReadTagWithCutoff(127);
    tag = p.first;
    if (!p.second) goto handle_unusual;
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // required int32 car_id = 1;
      case 1: {
        if (tag == 8) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &car_id_)));
          set_has_car_id();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(16)) goto parse_timestamp;
        break;
      }

      // required uint32 timestamp = 2;
      case 2: {
        if (tag == 16) {
         parse_timestamp:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::uint32, ::google::protobuf::internal::WireFormatLite::TYPE_UINT32>(
                 input, &timestamp_)));
          set_has_timestamp();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(24)) goto parse_lon;
        break;
      }

      // required int32 lon = 3;
      case 3: {
        if (tag == 24) {
         parse_lon:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &lon_)));
          set_has_lon();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(32)) goto parse_lat;
        break;
      }

      // required int32 lat = 4;
      case 4: {
        if (tag == 32) {
         parse_lat:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &lat_)));
          set_has_lat();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(40)) goto parse_head;
        break;
      }

      // required int32 head = 5;
      case 5: {
        if (tag == 40) {
         parse_head:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &head_)));
          set_has_head();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(48)) goto parse_speed;
        break;
      }

      // required int32 speed = 6;
      case 6: {
        if (tag == 48) {
         parse_speed:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, &speed_)));
          set_has_speed();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(61)) goto parse_x;
        break;
      }

      // required float x = 7;
      case 7: {
        if (tag == 61) {
         parse_x:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   float, ::google::protobuf::internal::WireFormatLite::TYPE_FLOAT>(
                 input, &x_)));
          set_has_x();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectTag(69)) goto parse_y;
        break;
      }

      // required float y = 8;
      case 8: {
        if (tag == 69) {
         parse_y:
          DO_((::google::protobuf::internal::WireFormatLite::ReadPrimitive<
                   float, ::google::protobuf::internal::WireFormatLite::TYPE_FLOAT>(
                 input, &y_)));
          set_has_y();
        } else {
          goto handle_unusual;
        }
        if (input->ExpectAtEnd()) goto success;
        break;
      }

      default: {
      handle_unusual:
        if (tag == 0 ||
            ::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          goto success;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, mutable_unknown_fields()));
        break;
      }
    }
  }
success:
  // @@protoc_insertion_point(parse_success:TrajPoint)
  return true;
failure:
  // @@protoc_insertion_point(parse_failure:TrajPoint)
  return false;
#undef DO_
}

void TrajPoint::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // @@protoc_insertion_point(serialize_start:TrajPoint)
  // required int32 car_id = 1;
  if (has_car_id()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(1, this->car_id(), output);
  }

  // required uint32 timestamp = 2;
  if (has_timestamp()) {
    ::google::protobuf::internal::WireFormatLite::WriteUInt32(2, this->timestamp(), output);
  }

  // required int32 lon = 3;
  if (has_lon()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(3, this->lon(), output);
  }

  // required int32 lat = 4;
  if (has_lat()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(4, this->lat(), output);
  }

  // required int32 head = 5;
  if (has_head()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(5, this->head(), output);
  }

  // required int32 speed = 6;
  if (has_speed()) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(6, this->speed(), output);
  }

  // required float x = 7;
  if (has_x()) {
    ::google::protobuf::internal::WireFormatLite::WriteFloat(7, this->x(), output);
  }

  // required float y = 8;
  if (has_y()) {
    ::google::protobuf::internal::WireFormatLite::WriteFloat(8, this->y(), output);
  }

  if (!unknown_fields().empty()) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        unknown_fields(), output);
  }
  // @@protoc_insertion_point(serialize_end:TrajPoint)
}

::google::protobuf::uint8* TrajPoint::SerializeWithCachedSizesToArray(
    ::google::protobuf::uint8* target) const {
  // @@protoc_insertion_point(serialize_to_array_start:TrajPoint)
  // required int32 car_id = 1;
  if (has_car_id()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(1, this->car_id(), target);
  }

  // required uint32 timestamp = 2;
  if (has_timestamp()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteUInt32ToArray(2, this->timestamp(), target);
  }

  // required int32 lon = 3;
  if (has_lon()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(3, this->lon(), target);
  }

  // required int32 lat = 4;
  if (has_lat()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(4, this->lat(), target);
  }

  // required int32 head = 5;
  if (has_head()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(5, this->head(), target);
  }

  // required int32 speed = 6;
  if (has_speed()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteInt32ToArray(6, this->speed(), target);
  }

  // required float x = 7;
  if (has_x()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteFloatToArray(7, this->x(), target);
  }

  // required float y = 8;
  if (has_y()) {
    target = ::google::protobuf::internal::WireFormatLite::WriteFloatToArray(8, this->y(), target);
  }

  if (!unknown_fields().empty()) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        unknown_fields(), target);
  }
  // @@protoc_insertion_point(serialize_to_array_end:TrajPoint)
  return target;
}

int TrajPoint::ByteSize() const {
  int total_size = 0;

  if (_has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    // required int32 car_id = 1;
    if (has_car_id()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int32Size(
          this->car_id());
    }

    // required uint32 timestamp = 2;
    if (has_timestamp()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::UInt32Size(
          this->timestamp());
    }

    // required int32 lon = 3;
    if (has_lon()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int32Size(
          this->lon());
    }

    // required int32 lat = 4;
    if (has_lat()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int32Size(
          this->lat());
    }

    // required int32 head = 5;
    if (has_head()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int32Size(
          this->head());
    }

    // required int32 speed = 6;
    if (has_speed()) {
      total_size += 1 +
        ::google::protobuf::internal::WireFormatLite::Int32Size(
          this->speed());
    }

    // required float x = 7;
    if (has_x()) {
      total_size += 1 + 4;
    }

    // required float y = 8;
    if (has_y()) {
      total_size += 1 + 4;
    }

  }
  if (!unknown_fields().empty()) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        unknown_fields());
  }
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void TrajPoint::MergeFrom(const ::google::protobuf::Message& from) {
  GOOGLE_CHECK_NE(&from, this);
  const TrajPoint* source =
    ::google::protobuf::internal::dynamic_cast_if_available<const TrajPoint*>(
      &from);
  if (source == NULL) {
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
    MergeFrom(*source);
  }
}

void TrajPoint::MergeFrom(const TrajPoint& from) {
  GOOGLE_CHECK_NE(&from, this);
  if (from._has_bits_[0 / 32] & (0xffu << (0 % 32))) {
    if (from.has_car_id()) {
      set_car_id(from.car_id());
    }
    if (from.has_timestamp()) {
      set_timestamp(from.timestamp());
    }
    if (from.has_lon()) {
      set_lon(from.lon());
    }
    if (from.has_lat()) {
      set_lat(from.lat());
    }
    if (from.has_head()) {
      set_head(from.head());
    }
    if (from.has_speed()) {
      set_speed(from.speed());
    }
    if (from.has_x()) {
      set_x(from.x());
    }
    if (from.has_y()) {
      set_y(from.y());
    }
  }
  mutable_unknown_fields()->MergeFrom(from.unknown_fields());
}

void TrajPoint::CopyFrom(const ::google::protobuf::Message& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void TrajPoint::CopyFrom(const TrajPoint& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool TrajPoint::IsInitialized() const {
  if ((_has_bits_[0] & 0x000000ff) != 0x000000ff) return false;

  return true;
}

void TrajPoint::Swap(TrajPoint* other) {
  if (other != this) {
    std::swap(car_id_, other->car_id_);
    std::swap(timestamp_, other->timestamp_);
    std::swap(lon_, other->lon_);
    std::swap(lat_, other->lat_);
    std::swap(head_, other->head_);
    std::swap(speed_, other->speed_);
    std::swap(x_, other->x_);
    std::swap(y_, other->y_);
    std::swap(_has_bits_[0], other->_has_bits_[0]);
    _unknown_fields_.Swap(&other->_unknown_fields_);
    std::swap(_cached_size_, other->_cached_size_);
  }
}

::google::protobuf::Metadata TrajPoint::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = TrajPoint_descriptor_;
  metadata.reflection = TrajPoint_reflection_;
  return metadata;
}


// @@protoc_insertion_point(namespace_scope)

// @@protoc_insertion_point(global_scope)