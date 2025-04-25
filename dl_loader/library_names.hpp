#pragma once

#ifndef DL_LOADER_LIBRARY_NAMES_HPP
#define DL_LOADER_LIBRARY_NAMES_HPP
#include <string>

#ifdef __APPLE__
static const std::string LIBRARY_NAME = "toupdate";
#endif

#ifdef __linux__
static const std::string LIBRARY_NAME = "toupdate";
#endif

#ifdef WIN32
static const std::string LIBRARY_NAME = "to_update";
#endif

#endif  // DL_LOADER_LIBRARY_NAMES_HPP