#pragma once
#ifndef IDL_LOADER_HPP
#define IDL_LOADER_HPP

#include <memory>

namespace dlloader {

template <typename T>
class IDLLoader {
 public:
  virtual ~IDLLoader() = default;

  virtual void DLOpenLib() = 0;
  virtual std::shared_ptr<T> DLGetInstance() = 0;
  virtual void DLCloseLib() = 0;
};
}  // namespace dlloader

#endif  // IDL_LOADER_HPP