import java.io.*;

import java.util.*;

import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.RecursiveAction;

import java.lang.Math;



// an RGB triple
class RGB {
    public int R, G, B;

    RGB(int r, int g, int b) {
	R = r;
	G = g;
	B = b;
    }

    public String toString() { return "(" + R + "," + G + "," + B + ")"; }

}

// code for creating a Gaussian filter
class Gaussian {

    protected static double gaussian(int x, int mu, double sigma) {
	return Math.exp( -(Math.pow((x-mu)/sigma,2.0))/2.0 );
    }

    public static double[][] gaussianFilter(int radius, double sigma) {
	int length = 2 * radius + 1;
	double[] hkernel = new double[length];
	for(int i=0; i < length; i++)
	    hkernel[i] = gaussian(i, radius, sigma);
	double[][] kernel2d = new double[length][length];
	double kernelsum = 0.0;
	for(int i=0; i < length; i++) {
	    for(int j=0; j < length; j++) {
		double elem = hkernel[i] * hkernel[j];
		kernelsum += elem;
		kernel2d[i][j] = elem;
	    }
	}
	for(int i=0; i < length; i++) {
	    for(int j=0; j < length; j++)
		kernel2d[i][j] /= kernelsum;
	}
	return kernel2d;
    }
}

// an object representing a single PPM image
class PPMImage {
    protected int width, height, maxColorVal;
    protected RGB[] pixels;

    PPMImage(int w, int h, int m, RGB[] p) {
	width = w;
	height = h;
	maxColorVal = m;
	pixels = p;
    }

	// parse a PPM file to produce a PPMImage
    public static PPMImage fromFile(String fname) throws FileNotFoundException, IOException {
	FileInputStream is = new FileInputStream(fname);
	BufferedReader br = new BufferedReader(new InputStreamReader(is));
	br.readLine(); // read the P6
	String[] dims = br.readLine().split(" "); // read width and height
	int width = Integer.parseInt(dims[0]);
	int height = Integer.parseInt(dims[1]);
	int max = Integer.parseInt(br.readLine()); // read max color value
	br.close();

	is = new FileInputStream(fname);
	    // skip the first three lines
	int newlines = 0;
	while (newlines < 3) {
	    int b = is.read();
	    if (b == 10)
		newlines++;
	}

	int MASK = 0xff;
	int numpixels = width * height;
	byte[] bytes = new byte[numpixels * 3];
        is.read(bytes);
	RGB[] pixels = new RGB[numpixels];
	for (int i = 0; i < numpixels; i++) {
	    int offset = i * 3;
	    pixels[i] = new RGB(bytes[offset] & MASK, bytes[offset+1] & MASK, bytes[offset+2] & MASK);
	}

	return new PPMImage(width, height, max, pixels);
    }

	// write a PPMImage object to a file
    public void toFile(String fname) throws IOException {
	FileOutputStream os = new FileOutputStream(fname);

	String header = "P6\n" + width + " " + height + "\n" + maxColorVal + "\n";
	os.write(header.getBytes());

	int numpixels = width * height;
	byte[] bytes = new byte[numpixels * 3];
	int i = 0;
	for (RGB rgb : pixels) {
	    bytes[i] = (byte) rgb.R;
	    bytes[i+1] = (byte) rgb.G;
	    bytes[i+2] = (byte) rgb.B;
	    i += 3;
	}
	os.write(bytes);
    }

    public PPMImage negate() 
	{
		Negate n = new Negate(width, height, maxColorVal, pixels);
		n.neg();
		return new PPMImage(n.width, n.height, n.maxColorVal, n.pixels);
        }

    public PPMImage mirrorImage() 
	{
		Mirror m = new Mirror(width, height, maxColorVal, pixels);
		m.mir();
		return new PPMImage(m.width, m.height, m.maxColorVal, m.pixels);
	}

    public PPMImage gaussianBlur(int radius, double sigma) 
    {
		GaussianBlur g = new GaussianBlur(width, height, maxColorVal, pixels, radius, sigma);
		g.gblur();
		return new PPMImage(g.width, g.height, g.maxColorVal, g.blurred);
    }
}

class Negate
{
	protected int width, height, maxColorVal;
	protected RGB pixels[];
	
	Negate(int w, int h, int m, RGB[] p) {
		width = w;
		height = h;
		maxColorVal = m;
		pixels = p;
	}

	public void neg() {
		for (RGB rgb : pixels) {
			rgb.R = maxColorVal - rgb.R;
			rgb.G = maxColorVal - rgb.G;
			rgb.B = maxColorVal - rgb.B;
		}
	}
}

class Mirror
{
	protected int width, height, maxColorVal;
	protected RGB pixels[];	

	Mirror(int w, int h, int m, RGB[] p) {
		width = w;
		height = h;
		maxColorVal = m;
		pixels = p;
	}

	public void mir() {
		RGB[] temp = new RGB[width];
		List<RGB> list = new ArrayList<RGB>();

		for (int i = 0; i < pixels.length; i+= width) {
			for (int j = 0 ; j < width; j++) {
				temp[j] = pixels[i+j];
			}
			for (int k = 0; k < width/2; k++) {
				RGB t = temp[k];
				temp[k] = temp[width-1-k];
				temp[width-1-k] = t;
			}
			for (int l = 0; l < width; l++) {
				list.add(temp[l]);
			}
		}
		pixels = list.toArray(new RGB[0]);
	}
}

class GaussianBlur 
{
	protected int width, height, maxColorVal, radius;
	protected RGB pixels[], blurred[];
	protected double[][] gfilter;
	
	GaussianBlur(int w, int h, int max, RGB[] p, int r, double sig) 
	{
		width = w;
		height = h;
		maxColorVal = max;
		pixels = p;
		radius = r;
		blurred = p;
		gfilter = Gaussian.gaussianFilter(r, sig);
	}

	public void gblur() 
	{
		for (int i = radius; i < height - radius; i++)
		{
			for (int j = radius; j < width - radius; j++)
			{
				double blurredR = 0.0;
				double blurredG = 0.0;
				double blurredB = 0.0;
				for (int di = -radius; di < radius + 1; di++)
				{
					for (int dj = -radius; dj < radius + 1; dj++)
					{
						blurredR += pixels[(i+di)*width+j+dj].R * gfilter[di+radius][dj+radius];
						blurredG += pixels[(i+di)*width+j+dj].G * gfilter[di+radius][dj+radius];
						blurredB += pixels[(i+di)*width+j+dj].B * gfilter[di+radius][dj+radius];
					}
				}
				blurred[i*width+j] = new RGB((int) blurredR, (int) blurredG, (int) blurredB);

			}
		}
	}
}

