#include "../include/Layer.h"
#include <shapelens/utils/FFT.h>

using namespace skylens;

ConvolutionLayer::ConvolutionLayer(double FOV, double pixsize, const PSF& psf) :
  psf(psf), pixsize(pixsize), L(int(floor(FOV/pixsize))),
  // automatically creates a single instance of LayerStack
  ls(SingleLayerStack::getInstance())
{
  Layer::z = 0;
  Layer::transparent = false;
  PAD = std::max(12,int(floor(0.01*L)));
  me = ls.insert(std::pair<double,Layer*>(z,this));

  // superimage does not exist yet:
  // we have to query all pixels of the image for values behind this Layer
  im.resize((L+2*PAD)*(L+2*PAD));
  im.grid = shapelens::Grid(-PAD,-PAD,L+2*PAD,L+2*PAD);
  im.grid.apply(shapelens::ScalarTransformation<double>(pixsize));

  LayerStack::iterator iter;
  for (unsigned long i = 0; i < im.size(); i++) {
    shapelens::Point<double> P = im.grid(i);
    im(i) = 0;
    iter = me;
    iter++;
    for (iter; iter != ls.end(); iter++) {
      if (iter->second->getRedshift() > 0) {
	std::string type = iter->second->getType();
	im(i) += iter->second->getFlux(P);
	if (type[0] == 'T')
	  break;
      }
    }
  }

  // now: convolve im with psf
  shapelens::ShapeletObject kernel = psf.getShape();
  shapelens::Object kernelObject = kernel.getModel();
  kernelObject.computeFluxCentroid();
  convolveImage(im,kernelObject);
}

double ConvolutionLayer::getFlux(const shapelens::Point<double>& P) const {
  if (!transparent) {
    return im.get(P);
  } else {
    double flux = 0;
    LayerStack::iterator iter = me;
    iter++; // next layer
    for (iter; iter != ls.end(); iter++) {
      if (iter->second->getRedshift() > 0) {
	std::string type = iter->second->getType();
	flux += iter->second->getFlux(P);
	if (type[0] == 'T')
	  break;
      }
    }
    return flux;
  }
}

std::string ConvolutionLayer::getType() const {
  return "TC";
}
// apply convolution to im
// therefore mask edges to avoid aliasing
void  ConvolutionLayer::convolveImage(shapelens::Image<double>& im, shapelens::Object& kernel) {
  for (int l=0; l < PAD; l++) {
    int j;
    double w = sin(l*M_PI_2/PAD);
    // low border
    int i = l;
    for (j=0; j<im.getSize(1); j++)
      im(i,j) *= w;
    // top border
    i = im.getSize(0) - l - 1;
    for (j=0; j<im.getSize(1); j++)
      im(i,j) *= w;
    // left border
    j = l;
    for (i=0; i<im.getSize(0); i++)
      im(i,j) *= w;
    // right border
    j = im.getSize(1) - l - 1;
    for (i=0; i<im.getSize(0); i++)
      im(i,j) *= w;
  }
  // compute FFT of im
  shapelens::FourierTransform2D fourier;
  shapelens::FFT::transform(im,fourier);
  
  // check size of kernel
  int n = im.getSize(0);
  int m = im.getSize(1);
  if (n != kernel.getSize(0)) {
    int N1 = kernel.getSize(0);
    int M1 = kernel.getSize(1); 
    if (N1%2==1) {
      N1++;
      kernel.resize_clear(N1*M1);
    }
    if (M1%2==1) {
      M1++;
      kernel.resize_clear(N1*M1);
    }
    int xmin = (N1-n)/2, xmax = N1 + (n-N1)/2, ymin = (M1-m)/2, ymax = M1+ (m-M1)/2;
    shapelens::Image<double> sub(xmax-xmin,ymax-ymin);
    kernel.slice(sub,shapelens::Point<int>(xmin,ymin),shapelens::Point<int>(xmax,ymax));
    kernel = sub;
    // compute new kernel::fourier
    kernel.computeFFT();
  }
  else if(kernel.fourier.getRealSize(0) == 0)
    kernel.computeFFT();

  // actual convolution
  shapelens::FFT::conv_multiply(fourier,kernel.fourier,fourier);
  // Transform back to real space and reorder:
  shapelens::FFT::transform(fourier,im);
  shapelens::FFT::reorder(im);
}
