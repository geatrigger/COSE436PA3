
#pragma once

#include <iostream>
#include "Viewer.h"

using namespace std;

Viewer OpenGL_Viewer;

void Initialize(void)
{
	OpenGL_Viewer.Initialize();
}


void Render(void)
{
	OpenGL_Viewer.Render();
}

void Reshape(int w, int h)
{
	OpenGL_Viewer.Reshape(w, h);
}

void Mouse(int mouse_event, int state, int x, int y)
{
	OpenGL_Viewer.Mouse(mouse_event, state, x, y);
}

void Motion(int x, int y)
{
	OpenGL_Viewer.Motion(x, y);
}

void PassiveMotion(int x, int y)
{
	OpenGL_Viewer.PassiveMotion(x, y);
}

void Keyboard(unsigned char key, int x, int y)
{
	OpenGL_Viewer.Keyboard(key, x, y);
}

void Update(int value)
{
	OpenGL_Viewer.Update();

	glutPostRedisplay();
	glutTimerFunc(10, Update, 0);
}

int main(int argc, char ** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
	glEnable(GL_DEPTH_TEST);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glutInitWindowPosition(400, 100);
	glutInitWindowSize(1024, 800);
	glutCreateWindow("Cloth Spring Simulation ");

	Initialize();
	glutTimerFunc(10, Update, 0);
	glutDisplayFunc(Render);

	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);
	glutPassiveMotionFunc(PassiveMotion);
	glutKeyboardFunc(Keyboard);
	glutReshapeFunc(Reshape);
	
	glutMainLoop();

	return 0;
}

