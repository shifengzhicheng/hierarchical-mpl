# 定义编译器和编译选项
CXX = g++
CXXFLAGS = -std=c++17 -fPIC -O2
LDFLAGS = -shared

# 查找Tcl库的包含路径和库路径
TCL_INCLUDE_PATH = /usr/include/tcl8.6
TCL_LIBRARY_PATH = /usr/lib/tcl8.6

# 定义源文件和目标文件
SRC_DIR = src
BUILD_DIR = build
LIB_DIR = $(BUILD_DIR)/lib
SRC_FILE = $(SRC_DIR)/graph_extension.cpp
TARGET = $(LIB_DIR)/libgraphextension.so

# 默认目标
all: $(TARGET)

# 生成共享库
$(TARGET): $(SRC_FILE)
	@mkdir -p $(LIB_DIR)
	$(CXX) $(CXXFLAGS) -I$(TCL_INCLUDE_PATH) $(LDFLAGS) -L$(TCL_LIBRARY_PATH) -ltcl -o $@ $<

# 清理生成的文件
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean