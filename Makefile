CXX = em++
CXXFLAGS = -O3 -std=c++20 -Wall -s WASM=1 -s EXPORT_ALL=1 --bind
TARGET=web/wasm.js

# Source files
SOURCES = code.cpp

# Build target
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

# Clean build
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: clean
