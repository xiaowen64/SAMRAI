#include	<stdio.h>
#include	<stdlib.h>
#include	<X11/Xlib.h>
#include	<GL/glx.h>
#include	<GL/gl.h>

/* Prototypes */

int main(int argc, char **argv);
void fatalError(char *);

int main(int argc, char **argv)
{
   Display		*dpy;
   XVisualInfo		*vi = NULL;

   int	sngBuf[] = {	GLX_RGBA,
				GLX_RED_SIZE, 1,
				GLX_GREEN_SIZE, 1,
				GLX_BLUE_SIZE, 1,
				GLX_DEPTH_SIZE, 12,
				None };

   int	dblBuf[] = {	GLX_RGBA,
				GLX_RED_SIZE, 1,
				GLX_GREEN_SIZE, 1,
				GLX_BLUE_SIZE, 1,
				GLX_DEPTH_SIZE, 12,
				GLX_DOUBLEBUFFER,
				None };

   Window		win;
   Bool		        doubleBuffer = True;
   int			dummy;

   if(!(dpy = XOpenDisplay(NULL)))
      fatalError("could not open X11 display");

   if(!glXQueryExtension(dpy, &dummy, &dummy))
      fatalError("X server has no OpenGL GLX extension\n");

   if(!(vi = glXChooseVisual(dpy, DefaultScreen(dpy), dblBuf))) {
      if(!(vi = glXChooseVisual(dpy, DefaultScreen(dpy), sngBuf)))
	 fatalError("no RGB visual with depth buffer\n");
      doubleBuffer = False;
   }
   if(vi->class != TrueColor)
      fatalError("TrueColor visual required for this program\n");

   return 0;
}

void fatalError(char *estr)
{
   fprintf(stderr,estr);
   exit(1);
}
