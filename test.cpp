
#include <iostream>
#include <string>
#include <SDL2/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

//#include <algebra_OO.h>
#include "algebra.hpp"

vec3_t box[8] = {{ 1.0, 1.0, 1.0}, {-1.0, 1.0, 1.0}, {-1.0,-1.0, 1.0}, {1.0,-1.0, 1.0},
                 { 1.0, 1.0,-1.0}, {-1.0, 1.0,-1.0}, {-1.0,-1.0, 1.0}, {1.0,-1.0,-1.0}};

class mainApp_c {
	private:
		int wnd_height, wnd_width;
		int wnd_posx, wnd_posy;
		
		SDL_Window *window;
		
		bool is_run;
		
	public:
		~mainApp_c();
		mainApp_c();

		void init_app();
		void looper();
};

mainApp_c::mainApp_c() {
	//std::cout << "HINT: Constructor call!" << std::endl;
	wnd_height = 640;
	wnd_width = 480;

	wnd_posx = SDL_WINDOWPOS_CENTERED;
	wnd_posy = SDL_WINDOWPOS_CENTERED;
	
	is_run = true;
}

mainApp_c::~mainApp_c() {

}

void mainApp_c::init_app() {
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		std::cout << "Error: Unable to init SDL; " << SDL_GetError() << std::endl;
		exit(1);
	}

	window = SDL_CreateWindow("Cube", wnd_posx, wnd_posy, 
							wnd_height, wnd_width, 
							SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);

	if(window == NULL) {
		std::cout << "Error: Unable to creaate window!" << std::endl;
		exit(1);
	}
	
	SDL_GLContext glcontext = SDL_GL_CreateContext(window); 

	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
	SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 6);
	SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 5);

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f); 
	glClearDepth(1.0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST); 
	glShadeModel(GL_SMOOTH);
}

void mainApp_c::looper() {
	while(is_run) {
		SDL_Event event; 
		
		while (SDL_PollEvent(&event)) {
			switch(event.type){ 
				case SDL_QUIT: 
					is_run = false;
					break;

				case SDL_KEYDOWN: 
					switch(event.key.keysym.sym) {
						case SDLK_ESCAPE: 
						is_run = false; 
						break;
					}
				break;
			}
		}
		
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(45.0f, (float) wnd_height / (float) wnd_width, 0.1f, 100.0f); 
		glTranslatef(0.0f, 0.0f, -8.0f);
		
		glMatrixMode(GL_MODELVIEW); 
		glLoadIdentity();
		
        glBegin(GL_QUADS);		// Рисуем куб

	    glColor3f(0.0f, 1.0f, 0.0f);
	    glVertex3fv(&box[0][0]);
	    glVertex3fv(&box[1][0]);
	    glVertex3fv(&box[2][0]);
	    glVertex3fv(&box[3][0]);
        
        glColor3f(1.0f, 1.0f, 0.0f);
	    glVertex3fv(&box[4][0]);
	    glVertex3fv(&box[5][0]);
	    glVertex3fv(&box[6][0]);
	    glVertex3fv(&box[7][0]);

	    glColor3f(0.0f, 1.0f, 1.0f);
	    glVertex3fv(&box[0][0]);
	    glVertex3fv(&box[1][0]);
	    glVertex3fv(&box[5][0]);
	    glVertex3fv(&box[4][0]);

	    glColor3f(1.0f, 0.0f, 0.0f);
	    glVertex3fv(&box[2][0]);
	    glVertex3fv(&box[3][0]);
	    glVertex3fv(&box[7][0]);
	    glVertex3fv(&box[6][0]);
	
	    glColor3f(0.0f, 0.0f, 1.0f);
	    glVertex3fv(&box[0][0]);
	    glVertex3fv(&box[4][0]);
	    glVertex3fv(&box[7][0]);
	    glVertex3fv(&box[3][0]);
	
	    glColor3f(1.0f, 0.0f, 1.0f);
	    glVertex3fv(&box[1][0]);
	    glVertex3fv(&box[5][0]);
	    glVertex3fv(&box[6][0]);
	    glVertex3fv(&box[2][0]);

	    glEnd();

		glFlush();
		SDL_GL_SwapWindow(window);	
	}
	
	SDL_Quit(); 
}

int main(int argc, char *argv[]) {
    std::cout << "ALGEBRA_TEST" << std::endl;
    mainApp_c app;

	app.init_app();

	app.looper();
}
