CXX := g++
CXXFLAGS := -std=c++14 -O3 -I. -MMD -MP -g 
LDFLAGS := 
SRCS := $(shell find . -name '*.cpp')
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)
TARGET := router


FLUTE_FILES := flute.h flute_int.h flute_malloc.h  
DATA_FILES := POWV9.dat POST9.dat

all: $(TARGET) $(FLUTE_FILES) $(DATA_FILES)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^


%.o: %.cpp $(FLUTE_FILES)
	$(CXX) $(CXXFLAGS) -Wno-sign-compare -Wno-unused-result -c $< -o $@



.PHONY: clean
clean:
	$(RM) $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)