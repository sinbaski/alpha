CXX := g++

ifdef DEBUG
CFLAGS := -Iinclude -c -g3 -fopenmp -Wall
OBJ_DIR := obj_debug
else
CFLAGS := -Iinclude -c -O3 -fopenmp -Wall
OBJ_DIR = obj
endif

CXXFLAGS := -std=c++11 $(CFLAGS)

LDFLAGS := -L/usr/lib/x86_64-linux-gnu -fopenmp
LDLIBS := -lgsl -lpthread -lm

CPPSRC := joule_heating.cpp

OBJ := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(CPPSRC))

.PHONY: prepare
prepare:
	@if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi


$(OBJ_DIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

jh.exe: $(OBJ)
	$(CXX) $^ $(LDLIBS) $(LDFLAGS) -o $@

.PHONY: deploy
deploy:
	mv -f jh.exe ..

clean:
	rm -f $(OBJ_DIR)/*.o jh.exe
