# Makefile
# Fluids Simulation
# N. Dommanget ndommanget@gmail.com



MACHINE = $(shell uname -s)

ifeq ($(MACHINE),Darwin)
	OPENGL_INC = -FOpenGL
	OPENGL_LIB = -framework OpenGL
else
	OPENGL_INC = -I/usr/X11R6/include
	OPENGL_LIB = -I/usr/lib64 -lGL -lGLU
endif

CXX=g++
COMPILER_FLAGS_WARN=-g -Wall 

INCLUDE= $(OPENGL_INC) -I../api -I.

SFML_LIB = -lsfml-graphics -lsfml-window -lsfml-system
LIBS= $(SFML_LIB) $(OPENGL_LIB)

OBJS = Object.o Builders.o Camera.o Simulation.o Scene.o Shaders.o Application.o Tools.o main.o glew.o

EXE=Demo


all : $(EXE) 

$(EXE) : $(OBJS)
	$(CXX) $(COMPILER_FLAGS_WARN) $^ $(LIBS) -o $@ 

-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CXX) -c $(COMPILER_FLAGS_WARN) $(INCLUDE) $<
	@$(CXX)  -MM $(INCLUDE) $< > $*.d

glew.o: glew/glew.c glew/glew.h glew/wglew.h glew/glxew.h
	$(CXX) -DGLEW_NO_GLU $(COMPILER_FLAGS) -o $@ -c $<

.PHONY : clean bigclean
clean::
	rm -f *.o *~ *.d
bigclean : clean
	rm -f $(EXE)

