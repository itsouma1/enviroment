# Build configuration for the climate analysis tool.
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
LDFLAGS  ?=
BIN      := climate

SRC      := $(wildcard src/*.cpp)
OBJ      := $(SRC:.cpp=.o)
OBJ_NOMAIN := $(filter-out src/main.o,$(OBJ))

all: $(BIN) tools/generate_demo

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJ)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

tools/generate_demo: tools/generate_demo.cpp $(OBJ_NOMAIN)
	$(CXX) $(CXXFLAGS) -Isrc -o $@ tools/generate_demo.cpp $(OBJ_NOMAIN)

tests/test_all: tests/test_all.cpp $(OBJ_NOMAIN)
	$(CXX) $(CXXFLAGS) -Isrc -o $@ tests/test_all.cpp $(OBJ_NOMAIN)

test: tests/test_all
	./tests/test_all

data: tools/generate_demo
	./tools/generate_demo data

demo: $(BIN) data
	./$(BIN) query "is the arctic changing?" --data data

ui: $(BIN)
	./$(BIN) ui --port 8080

clean:
	rm -f $(OBJ) $(BIN) tools/generate_demo tests/test_all

.PHONY: all test data demo ui clean
