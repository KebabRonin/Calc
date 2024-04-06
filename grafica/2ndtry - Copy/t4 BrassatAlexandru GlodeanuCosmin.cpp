#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <algorithm>

#include "glut.h"
// the size of the window measured in pixels
#define dim 300
#define PI 3.14159265
unsigned char prevKey;
int line_width = 1;

struct Point { float x, y; };

class CartesianGrid {
public:
    CartesianGrid(unsigned int width, unsigned int height) : width(width), height(height) {
        r = 0.6 / ((width > height ? width : height) - 1);
    }
    void display() {
        glColor3f(1.0, 0.1, 0.1);
        glBegin(GL_LINES);
        for (int i = 0; i < height; ++i) {
            glVertex2f(-0.9 , 0.9 - (i * 1.8) / (height-1));
            glVertex2f(0.9, 0.9 - (i * 1.8) / (height - 1));
            for (int j = 0; j < width; ++j) {
                glVertex2f(0.9 - (j * 1.8) / (width - 1), -0.9);
                glVertex2f(0.9 - (j * 1.8) / (width - 1),  0.9);
            }
        }
        glEnd();
    }
    void writePixel(unsigned int i, unsigned int j) {
        //printf("Wrote to %d, %d\n", i, j);
        float x = -0.9 + (i * 1.8) / (width  - 1),
              y = -0.9 + (j * 1.8) / (height - 1);
        if (!(i < width && j < height)) {
            //printf("%d, %d iese din grid\n", i, j);
            //exit(0);
            return;
        }

        glColor3f(0.0, 0.0, 0.0); // negru
        glBegin(GL_POLYGON);
            for (float angle = 2 * PI; angle > 0 ; angle -= 2 * PI / 36)
            {
                glVertex2f(x + cos(angle) * r, y + sin(angle) * r);
            }
        glEnd();
    }
    Point PixelToScreenPoint(unsigned int i, unsigned int j) {
        Point p = { -0.9 + (i * 1.8) / (width - 1), -0.9 + (j * 1.8) / (height - 1) };
        return p;
    }
    Point ScreenPointToGridPoint(Point p) {
        Point pg = { (p.x - 0.1 + 1)/2.0 * (width+1), (p.y - 0.1 + 1)/2.0 * (height+1) }; // ????
        return pg;
    }
private:
    unsigned int width, height;
    float r;
};


void ScanConvertSegments3(Point p1, Point p2, CartesianGrid g, unsigned int grosime) {
    // the initial value of the decision variable
    // dx, dy are constants - see above
    if (p1.x > p2.x) {
        std::swap(p1, p2);
    }
    int sign = (p2.y - p1.y) / (p2.x - p1.x) < 0 ? -1 : 1;
    //printf("ddd %d %f %f, %f %f\n", sign, p1.x, p1.y, p2.x, p2.y);
    int dx = (p2.x - p1.x);
    int dy = (p2.y - p1.y);
    int d = 2 * dy - dx * sign;
    int dE = 2 * dy;
    int dNE = 2 * (dy - dx * sign);

    //SE: 2(a(x+2)+b(y-1/2)+c)=
    //    2(a(x+1)+b(y+1/2)+c)+2(a-b)
    float x = p1.x, y = p1.y;
    g.writePixel(x, y);
    for (int i = 1; i < grosime; ++i) {
        g.writePixel(x, y + i);
        g.writePixel(x, y - i);
    }
    while (x < p2.x)
    {
        //printf("11 %d %d %d %f %f\n", dE, dNE, d, x, y);
        if (d * sign <= 0) { /* select E */ 
            d += dE; 
            x++; 
        }
        else { /* select NE */ 
            d += dNE; 
            x++, y+= sign; 
        }
        //printf("22 %d %d %d %f %f\n", dE, dNE, d, x, y);
    //    M = M . (x, y);
        g.writePixel(x, y);
        for (int i = 1; i < grosime; ++i) {
            g.writePixel(x, y + i);
            g.writePixel(x, y - i);
        }
    }
}


void Display1() {
    CartesianGrid g = CartesianGrid(16, 16);
    g.display();
    Point p1, p2;
    glLineWidth(line_width+3);
    glColor3f(1.0f, 0.1f, 0.1f); // rosu
    glBegin(GL_LINES);
        p1 = g.PixelToScreenPoint(0, 0);
        p2 = g.PixelToScreenPoint(15, 7);
        //printf("%f %f, %f %f\n", p1.x, p1.y, p2.x, p2.y);

        glVertex2f(p1.x, p1.y);
        glVertex2f(p2.x, p2.y);
    glEnd();
    glLineWidth(line_width);
    ScanConvertSegments3(g.ScreenPointToGridPoint(p1), g.ScreenPointToGridPoint(p2), g, 1);

    glLineWidth(line_width+3);
    glColor3f(1.0f, 0.1f, 0.1f); // rosu
    glBegin(GL_LINES);
        p1 = g.PixelToScreenPoint(0, 15);
        p2 = g.PixelToScreenPoint(15, 10);
        //printf("%f %f, %f %f\n", p1.x, p1.y, p2.x, p2.y);

        glVertex2f(p1.x, p1.y);
        glVertex2f(p2.x, p2.y);
    glEnd();
    glLineWidth(line_width);
    ScanConvertSegments3(g.ScreenPointToGridPoint(p1), g.ScreenPointToGridPoint(p2), g, 2);
}

void ScanConvertCircle4(float r, CartesianGrid g, int grosime) {
    int x = 0, y = r;
    int d = 1 - r;
    int dE = 3, dSE = -2 * r + 5;
    //ScanConvertCircle3_Aux(x, y, M);
    // not used because only one octant
    g.writePixel(y, x);
    for (int i = 1; i < grosime; ++i) {
        g.writePixel(y + i, x);
        g.writePixel(y - i, x);
    }
    while (y > x) {
        if (d < 0)
        {
            d += dE;
            dE += 2;
            dSE += 2;
        }
        else
        {
            d += dSE;
            dE += 2;
            dSE += 4;
            y--;
        }
        x++;
        g.writePixel(y, x);
        for (int i = 1; i < grosime; ++i) {
            g.writePixel(y + i, x);
            g.writePixel(y - i, x);
        }
        //ScanConvertCircle3_Aux(x, y, M);
        // not used because only one octant
    }
}

void Display2() {
    CartesianGrid g = CartesianGrid(16, 16);
    g.display();
    Point p;
    glLineWidth(line_width+3);
    glColor3f(1.0f, 0.1f, 0.1f); // rosu
    glBegin(GL_LINE_LOOP);
        p = g.PixelToScreenPoint(0, 0);
        float r = g.PixelToScreenPoint(0, 13).y - p.y;
        for (float angle = 2 * PI; angle > 0; angle -= 2 * PI / 72)
        {
            glVertex2f(p.x + cos(angle) * r, p.y + sin(angle) * r);
        }
    glEnd();
    glLineWidth(line_width);

    ScanConvertCircle4(13, g, 2);
}

void Display(void) 
{
  switch(prevKey) 
  {
    case '0':
      glClear(GL_COLOR_BUFFER_BIT);
      break;
    case '1':
      glClear(GL_COLOR_BUFFER_BIT);
      Display1();
      break;
    case '2':
        glClear(GL_COLOR_BUFFER_BIT);
        Display2();
        break;
    default:
      break;
  }

  glFlush();
}

void Reshape(int w, int h) 
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

void KeyboardFunc(unsigned char key, int x, int y) 
{
   prevKey = key;
   if (key == 'q')
      exit(0);
   glutPostRedisplay();
}

void MouseFunc(int button, int state, int x, int y) 
{

}

void Init(void) {

    glClearColor(1.0, 1.0, 1.0, 1.0);

    glLineWidth(1);

    glPolygonMode(GL_FRONT, GL_LINE);
}

int main(int argc, char** argv) 
{
  glutInit(&argc, argv);

  glutInitWindowSize(dim, dim);

  glutInitWindowPosition(100, 100);

  glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);

  glutCreateWindow (argv[0]);

  Init();

  glutReshapeFunc(Reshape);

  glutKeyboardFunc(KeyboardFunc);

  glutMouseFunc(MouseFunc);

  glutDisplayFunc(Display);

  glutMainLoop();

  return 0;
}


