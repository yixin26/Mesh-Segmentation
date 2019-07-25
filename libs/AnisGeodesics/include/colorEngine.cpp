#include "colorEngine.h"


void ColorEngine::colorHSV(unsigned n, std::vector<color> &colors)
{
	unsigned binSize = (unsigned)(ceil((float)(n/6.0f)));

	// from red to yellow
	// increase green intensity
	for (unsigned int i=0; i< binSize; ++i)
	{
		color c;
		c.r = 1.0f;			
		c.g = ((float)i/(float)binSize);
		c.b = 0.0f;
		colors.push_back(c);
	}

	// from yellow to green
	// decrease red intensity
	for (unsigned int i=binSize; i >0 ; --i)
	{
		color c;
		c.r = ((float)i/(float)binSize);
		c.g = 1.0f;
		c.b = 0.0f;
		colors.push_back(c);
	}

	// from green to cyan
	// increase blue intensity
	for (unsigned int i=0; i< binSize; ++i)
	{
		color c;
		c.r = 1.0f;
		c.g = 0.0f;
		c.b =((float)i/(float)binSize);
		//color c(0.0f,1.0f, (i/binSize));
		colors.push_back(c);
	}

	// from cyan to blue
	// decrease green intensity
	for (unsigned int i=binSize; i >0 ; --i)
	{
		color c;
		c.r = 0.0f;
		c.g = ((float)i/(float)binSize);
		c.b = 1.0f;
		//color c(0.0f, (i/binSize), 1.0f);
		colors.push_back(c);
	}

	// from blue to magenta
	// increase red intensity
	for (unsigned int i=0; i< binSize; ++i)
	{
		color c;
		c.r = ((float)i/(float)binSize);
		c.g = 0.0f;
		c.b = 1.0f;
		//color c((i/binSize), 0.0f, 1.0f);
		colors.push_back(c);
	}

	// from magenta to red
	// decrease blue intensity
	for (unsigned int i=binSize; i > 0 ; --i)
	{
		color c;
		c.r = 1.0f;
		c.g = 0.0f;
		c.b = ((float)i/(float)binSize);
		//color c(1.0f, 0.0, (i/binSize));
		colors.push_back(c);
	}
}
void ColorEngine::colorHeat(unsigned n, std::vector<color> &colors)
{
	unsigned binSize = (unsigned)(ceil((float)(n/4.0f)));

	// from R to Y		
	for (unsigned int i=0; i< binSize; ++i)
	{
		color c;
		c.r = 1.0f;			
		c.g = ((float)i/(float)binSize);
		c.b = 0.0f;
		colors.push_back(c);
	}

	// from Y to Green
	for (unsigned int i=binSize; i > 0 ; --i)
	{
		color c;
		c.r = ((float)i/(float)binSize);
		c.g = 1.0f;
		c.b = 0.0f;
		colors.push_back(c);
	}

	// from Green to C
	for (unsigned int i=0; i< binSize; ++i)
	{
		color c;
		c.r = 0.0f;
		c.g = 1.0f;
		c.b = ((float)i/(float)binSize);
		colors.push_back(c);
	}

	// from C to B
	for (unsigned int i=binSize; i >0 ; --i)
	{
		color c;
		c.r = 0.0f;
		c.g = ((float)i/(float)binSize);
		c.b = 1.0f;
		colors.push_back(c);
	}
}
void ColorEngine::HslToRgb(std::vector<double>& mag, std::vector<color> &colors)
{
	HSL hsl; hsl.Y=1.0;hsl.Z=0.5;
	color tcolor;
	float& r=tcolor.r;
	float& g=tcolor.g;
	float& b=tcolor.b;

	for(unsigned i=0;i<mag.size();i++)
	{
		hsl.X=mag[i];

		if (hsl.Y == 0.0f)
			r = g = b = hsl.Z;
		else
		{
			double q = hsl.Z < 0.5f ? hsl.Z * (1.0f + hsl.Y) : hsl.Z + hsl.Y - hsl.Z * hsl.Y;
			double p = 2.0f * hsl.Z - q;
			r = HueToRgb(p, q, hsl.X + 1.0f / 3.0f);
			g = HueToRgb(p, q, hsl.X);
			b = HueToRgb(p, q, hsl.X - 1.0f / 3.0f);
		}
/*
		if(tcolor.r ==0 &&tcolor.g ==0 &&tcolor.b ==0 ) {
			cout<<"-"<<i<<":"<<mag[i]<<endl;
		}
*/
		colors.push_back(tcolor);
	}
}
double ColorEngine::HueToRgb(double p, double q, double t)
{
	if (t < 0.0) t += 1.0f;
	if (t > 1.0) t -= 1.0f;
	if (t < 1.0 / 6.0) return p + (q - p) * 6.0 * t;
	if (t < 1.0 / 2.0) return q;
	if (t < 2.0 / 3.0) return p + (q - p) * (2.0/3.0 - t) * 6.0;
	return p;
}