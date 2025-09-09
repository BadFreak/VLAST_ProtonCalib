cmake_minimum_required(VERSION 3.10)
project(VLAST_ProtonCalib)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 查找 ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

set(GENFIT_DIR "/mnt2/USTC/jxwang/Prefix/GenFit")

include_directories(${GENFIT_DIR}/include)
link_directories(${GENFIT_DIR}/lib64)

# 头文件目录
include_directories(include)

# 源文件
set(SRC_FILES src/main.cxx src/ProtonCalib.cxx)

# 生成可执行文件
add_executable(VLAST_ProtonCalib ${SRC_FILES})

# 链接 ROOT 库
target_link_libraries(VLAST_ProtonCalib genfit2 ${ROOT_LIBRARIES})
