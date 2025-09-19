#ifndef PRATIK_2_MYNODE__VISIBILITY_CONTROL_H_
#define PRATIK_2_MYNODE__VISIBILITY_CONTROL_H_
#if defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
    #define PRATIK_2_MYNODE_EXPORT __attribute__ ((dllexport))
    #define PRATIK_2_MYNODE_IMPORT __attribute__ ((dllimport))
  #else
    #define PRATIK_2_MYNODE_EXPORT __declspec(dllexport)
    #define PRATIK_2_MYNODE_IMPORT __declspec(dllimport)
  #endif
  #ifdef PRATIK_2_MYNODE_BUILDING_LIBRARY
    #define PRATIK_2_MYNODE_PUBLIC PRATIK_2_MYNODE_EXPORT
  #else
    #define PRATIK_2_MYNODE_PUBLIC PRATIK_2_MYNODE_IMPORT
  #endif
  #define PRATIK_2_MYNODE_PUBLIC_TYPE PRATIK_2_MYNODE_PUBLIC
  #define PRATIK_2_MYNODE_LOCAL
#else
  #define PRATIK_2_MYNODE_EXPORT __attribute__ ((visibility("default")))
  #define PRATIK_2_MYNODE_IMPORT
  #if __GNUC__ >= 4
    #define PRATIK_2_MYNODE_PUBLIC __attribute__ ((visibility("default")))
    #define PRATIK_2_MYNODE_LOCAL  __attribute__ ((visibility("hidden")))
  #else
    #define PRATIK_2_MYNODE_PUBLIC
    #define PRATIK_2_MYNODE_LOCAL
  #endif
  #define PRATIK_2_MYNODE_PUBLIC_TYPE
#endif
#endif  // PRATIK_2_MYNODE__VISIBILITY_CONTROL_H_
// Generated 17-Sep-2025 17:30:18
// Copyright 2019-2020 The MathWorks, Inc.
