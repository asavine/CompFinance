#pragma once
#ifndef IDL_LOADER_UNIX_HPP
#define IDL_LOADER_UNIX_HPP

#if defined(__linux__) || defined(__APPLE__)
#include <dlfcn.h>

#include <iostream>
#include <string>

#include "../IDLLoader.hpp"

namespace dlloader {

template <typename T>
class DLLoaderUnix : public IDLLoader<T> {
  void *handle;
  std::string pathToLib;
  std::string allocClassSymbol;
  std::string deallocClassSymbol;

 public:
  DLLoaderUnix(const std::string &pathToLib,
               const std::string &allocClassSymbol,
               const std::string &deallocClassSymbol)
      : handle(nullptr),
        pathToLib(pathToLib),
        allocClassSymbol(allocClassSymbol),
        deallocClassSymbol(deallocClassSymbol) {}

  ~DLLoaderUnix() = default;

  void DLOpenLib() override {
    handle = dlopen(pathToLib.c_str(), RTLD_NOW | RTLD_LAZY);
    if (!handle) {
      //   throw std::runtime_error(dlerror());
      std::cerr << "Error opening library: " << dlerror() << std::endl;
    }
  }

  std::shared_ptr<T> DLGetInstance() override {
    using allocClass = T *(*)();
    using deallocClass = void (*)(T *);

    auto allocFunc =
        reinterpret_cast<allocClass>(dlsym(handle, allocClassSymbol.c_str()));

    auto deleteFunc = reinterpret_cast<deallocClass>(
        dlsym(handle, deallocClassSymbol.c_str()));
    if (!allocFunc || !deleteFunc) {
      DLCloseLib();
      std::cerr << "Error loading symbols: " << dlerror() << std::endl;
    }

    return std::shared_ptr<T>(allocFunc(),
                              [deleteFunc](T *ptr) { deleteFunc(ptr); });
  }

  void DLCloseLib() override {
    if (dlclose(handle) != 0) {
      std::cerr << "Error closing library: " << dlerror() << std::endl;
    }
  }
};
}  // namespace dlloader
#endif  // defined(__linux__) || defined(__APPLE__)

#endif  // IDL_LOADER_UNIX_HPP