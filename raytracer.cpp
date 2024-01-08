#include <stdio.h>
#include <valarray>
/*---		Output Function		---*/
// Output in P6 format, a binary file containing:
// P6
// ncolumns nrows
// Max colour value
// colours in binary format thus unreadable
void save_imageP6(int Width, int Height, char* fname,unsigned char* pixels) {
  FILE *fp;
  const int maxVal=255; 
  
  printf("Saving image %s: %d x %d\n", fname,Width,Height);
  fp = fopen(fname,"wb");
  if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
  }
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", Width, Height);
  fprintf(fp, "%d\n", maxVal);

  for(int j = 0; j < Height; j++) {
		  fwrite(&pixels[j*Width*3], 3,Width,fp);
  }

  fclose(fp);
}

// Output in P3 format, a text file containing:
// P3
// ncolumns nrows
// Max colour value (for us, and usually 255)
// r1 g1 b1 r2 g2 b2 .....
void save_imageP3(int Width, int Height, char* fname,unsigned char* pixels) {
	FILE *fp;
	const int maxVal=255;
	
	printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"w");
	if (!fp) {
		printf("Unable to open file '%s'\n",fname);
		return;
	}
	fprintf(fp, "P3\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);
	
	int k = 0 ;
	for(int j = 0; j < Height; j++) {
		
		for( int i = 0 ; i < Width; i++)
		{
			fprintf(fp," %d %d %d", pixels[k],pixels[k+1],pixels[k+2]) ;
			k = k + 3 ;
		}
		fprintf(fp,"\n") ;
	}
	fclose(fp);
}

/* --- MAIN PROGRAM --- */
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
using namespace std;

double near[5];							// Near plane. near, left, right, bottom, top
int resolution[2];					// Resolution. Width x Length
double eye[3] = {0.0, 0.0, 0.0};				// Eye position
double backgroundColor[3];				// Background RGB values
double ambientIntensity[3];				// Intensity RGB values
double maxDepth = 100;
char output_name[21];					// Output file name	(21 characters no spaces)

class Ray {	// compute pixel in world coordinates, and also ray from eye -> pixel
public:
	double point[3];		// world coordinate of pixel, offset needed???
	double origin[3];		//added these two for raytracing
    double direction[3];

    Ray(const double origin[3], const double direction[3]) {	// for raytracing, point plus a vector
        for (int i = 0; i < 3; i++) {
            this->origin[i] = origin[i];
            this->direction[i] = direction[i];
			this->point[i] = direction[i] + origin[i];
        }
    }
	Ray(double u, double v) {	// u, v are pixel coordinates starting from top left.
		point[0] = near[2] * ((2.0*u)/(double)resolution[0]-1.0);
		point[1] = near[4] * ((2.0*v)/(double)resolution[1]-1.0);
		point[2] = -near[0];

		for (int i = 0; i < 3; i++) {
			origin[i] = 0;
			direction[i] = point[i];
		}
		normalize(direction);
	}
private:
	void normalize(double vec[3]) {
		double length = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
		if (length != 0.0) {
			for (int i = 0; i < 3; i++) {
				vec[i] /= length;
			}
		}
	}
};

class Light {
public:
	char name[21];
	double position[3];					// Position x,y,z
	double intensity[3];					// Light source R,G,B
	
	Light(const char* newName, double x, double y, double z, double r, double g, double b)
	{
		position[0] = (x); position[1] = (y); position[2] = (z);
		intensity[0] = (r); intensity[1] = (g); intensity[2] = (b);
		strncpy(name, newName, sizeof(name) -1);
		name[sizeof(name) - 1] = '\0';
	}

};

class Sphere {
public:
	char name[21];
	double position[3];					// Position X,Y,Z
	double scale[3];						// Scale X,Y,Z
	double color[3];						// Color R,G,B	
	double phong[5];			// 0 Ambience, 1 diffuse, 2 specular, 3 reflection, 4 specular exponent
	Sphere(const char* newName, double x, double y, double z, double i, double j, double k, double r, double g, double b, double a, double d, double s, double rf, double n)
		{
		position[0] = (x); position[1] = (y); position[2] = (z); scale[0] = (i); scale[1] = (j); scale[2] = (k);
		color[0] = (r); color[1] = (g); color[2] = (b);
		phong[0] = (a); phong[1] = (d); phong[2] = (s); phong[3] = (rf); phong[4] = (n);
		
		strncpy(name, newName, sizeof(name) -1);
		name[sizeof(name) - 1] = '\0';
	}
	bool intersect(Ray ray, double& t) const {
        // ray-sphere intersection formula. oc is the vector from the sphere origin to the ray origin (near plane)
        double oc[3] = {((ray.point[0] - position[0])/scale[0]),
			((ray.point[1] - position[1])/scale[1]),
			((ray.point[2] - position[2])/scale[2])};
		double rayDir[3] = {
			ray.point[0] /scale[0],
			ray.point[1] /scale[1],
			ray.point[2] /scale[2]
    	};
		double a = dot(rayDir, rayDir);
		double b = 2.0 * dot(oc, rayDir);
		double c = dot(oc, oc) - 1.0;

        // calculate discriminant
        double discriminant = (b * b) - (4 * a * c);

        if (discriminant < 0) {
            // no intersection
            return false;
        } else {
            // compute the two possible solutions for t
            double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
            double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
            // check if the intersection point is in front of the camera
            if (t1 > 0 || t2 > 0) {
                // choose the smaller positive solution
                t = (t1 < t2) ? t1 : t2;
                return true;
            } else {
                // both solutions are behind the camera
                return false;
            }
        }
    }
	bool intersect(const double startPoint[3], const double direction[3], double& t) const {
		// ray-sphere intersection formula. oc is the vector from the sphere origin to the ray origin
		double oc[3] = {(startPoint[0] - position[0]) / scale[0],
						(startPoint[1] - position[1]) / scale[1],
						(startPoint[2] - position[2]) / scale[2]};

		double rayDir[3] = {direction[0] / scale[0],
							direction[1] / scale[1],
							direction[2] / scale[2]};

		double a = dot(rayDir, rayDir);
		double b = 2.0 * dot(oc, rayDir);
		double c = dot(oc, oc) - 1.0;

		// calculate discriminant
		double discriminant = (b * b) - (4 * a * c);

		if (discriminant < 0) {
			// no intersection
			return false;
		} else {
			// compute the two possible solutions for t
			double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
			double t2 = (-b + sqrt(discriminant)) / (2.0 * a);

			// check if the intersection point is in front of the camera
			if (t1 > 0 || t2 > 0) {
				// choose the smaller positive solution
				t = (t1 < t2) ? t1 : t2;
				return true;
			} else {
				// both solutions are behind the camera
				return false;
			}
		}
	}
	double* pixelColor(vector<Light*> lights, double hitpoint[3], double pixelColor[3], vector<Sphere*> spheres){
		double normal[3];
		getNormal(hitpoint, normal);
		for (int i = 0; i < 3; i++) {	// ambient component
			pixelColor[i] = phong[0] * ambientIntensity[i] * color[i];
		}

		double directionToLight[3];	// diffuse component
		for (int j = 0; j < lights.size(); j++) {	// for each light
			for (int i = 0; i < 3; i++) {	// for each pixel channel
				for( int k = 0; k < 3; k++) {
					directionToLight[k] = lights.at(j)->position[k] - hitpoint[k];
				}
				normalize(directionToLight);
				double dotP = std::max(0.0, dot(normal, directionToLight));
				pixelColor[i] += (phong[1] * lights.at(j)->intensity[i] * dotP * color[i]);

				pixelColor[i] = std::min(pixelColor[i], 1.0);	// max out at 1
			}
		}

		
		double directionToEye[3]; // specular component
		double reflected[3]; // reflected ray
		for (int j = 0; j < lights.size(); j++) {	// for each light
			if(!isInShadow(*lights.at(j), hitpoint, spheres))
			{
				for (int i = 0; i < 3; i++) {	// for each pixel channel
					for( int k = 0; k < 3; k++) {	// calculate direction to eye vector
						directionToEye[k] = eye[k] - hitpoint[k];
						directionToLight[k] = lights.at(j)->position[k] - hitpoint[k];
					}

					normalize(directionToEye);
					normalize(directionToLight);
					getReflected(hitpoint, normal, directionToLight, reflected);

					double dotP = max(0.0, dot(reflected, directionToEye));	// should be reflected ray and direction to eye
					pixelColor[i] += (phong[2] * lights.at(j)->intensity[i] * pow(dotP, phong[4]));

					pixelColor[i] = std::min(pixelColor[i], 1.0);	// max out at 1
			
				}
			}
	
		}

		return pixelColor;
	}

	bool isInShadow(const Light light, const double intersectionPoint[3], vector<Sphere*> spheres) const {

        double shadowRayDirection[3];
		double tempIntersection[3];
        for (int i = 0; i < 3; ++i) {
            shadowRayDirection[i] = light.position[i] - intersectionPoint[i];
			tempIntersection[i] = intersectionPoint[i] + 0.000001;
        }

        double shadowRayLength = std::sqrt(shadowRayDirection[0] * shadowRayDirection[0] +
                                           shadowRayDirection[1] * shadowRayDirection[1] +
                                           shadowRayDirection[2] * shadowRayDirection[2]);
        for (int i = 0; i < 3; ++i) {	// normalize
            shadowRayDirection[i] /= shadowRayLength;
        }
        Ray shadowRay(tempIntersection, shadowRayDirection);

        // check for intersections with other objects in the scene
        for (const auto& otherSphere : spheres) {
            double t = maxDepth;
            if (otherSphere->intersect(shadowRay.origin, shadowRay.direction, t)) {
				return true;
            }
        }

        // not in shadow
        return false;
    }

private:
	double dot(const double a[3], const double b[3]) const {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
	void increment(double* point) {
		for (int i = 0; i < 3; i++) {
			point[i] += 0.000001;
		}

	}
	void getNormal(double intersectionPoint[3], double normal[3]) const {
        for (int i = 0; i < 3; ++i) {
            normal[i] = (intersectionPoint[i] - position[i]);
        }
        double length = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]); // normalize
        for (int i = 0; i < 3; ++i) {
            normal[i] /= length;
        }
		
		double directionToViewer[3];
		for (int i = 0; i < 3; ++i) {
			directionToViewer[i] = eye[i] - intersectionPoint[i];
		}
		
		double dotProduct = dot(normal, directionToViewer);	// make sure normal is in right direction
		if (dotProduct < 0.0) {
			// Flip the normal if it points inward
			for (int i = 0; i < 3; ++i) {
				normal[i] *= -1.0;
			}
		}
    }

	void getReflected(double intersectionPoint[3], double normal[3], double toLight[3], double reflected[3]) const {
		
		double length = sqrt(toLight[0] * toLight[0] + toLight[1] * toLight[1] + toLight[2] * toLight[2]); // normalize
        for (int i = 0; i < 3; ++i) {
            toLight[i] /= length;
        }
		double dotProduct = 2.0 * dot(toLight, normal);

		for (int i = 0; i < 3; i++) {
				reflected[i] = dotProduct * normal[i] - toLight[i];
		}
	}

	void normalize(double vec[3]) {
		double length = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
		if (length != 0.0) {
			for (int i = 0; i < 3; i++) {
				vec[i] /= length;
			}
		}
	}
};


/* Input Format...
 * NEAR n
 * LEFT l
 * RIGHT r
 * BOTTOM b
 * TOP t
 * RES x y
 * SPHERE name posx posy posz sclx scly sclz r g b ka kd ks kr n
 * LIGHT name posx posy posz ir ig ib
 * BACK r g b
 * AMBIENT ir ig ib
 * OUTPUT name
 */


int main(int argc, char *argv[]) {
	
	std::vector<Sphere*> spheres;
	std::vector<Light*> lights;

	/* SETUP CODE */
	if (argc != 2) {	// Check for proper usage
		std::cout << "Use: " << argv[0] << "input.txt" << std::endl;
		return 1;
	}

	std::ifstream file(argv[1]);

	if (!file.is_open()) {		//  check for file opening issues
		std::cout << "Error opening file " << argv[1] << "." << std::endl;
		return 1;
	}

	std::string command;
	int tmp = 0;

	while (file >> command) {	// Setup Code
		std::cout << "Command: " << command << std::endl;  // Debug print

		if (command == "NEAR" || command == "LEFT" || command == "RIGHT" || command == "BOTTOM" || command == "TOP") {
			file >> near[tmp];
			if (near[0] < 0) { near[0] = -1 * near[0]; }
			std::cout << "Near Value: " << near[tmp - 1] << std::endl;  // Debug print
			tmp++;
		}
		else if (command == "RES") {
			file >> resolution[0] >> resolution [1];
			std::cout << "Resolution: " << resolution[0] << " " << resolution[1] << std::endl;  // Debug print
		}
		else if (command == "SPHERE") {
			// buffer values to pass into constructor
			char name[21]; file >> name;
			double x, y, z; file >> x; file >> y; file >> z;
			double i, j, k; file >> i; file >> j; file >> k;
			double r, g, b; file >> r; file >> g; file >> b;
			double a, d, s, rf, n; file >> a; file >> d; file >> s; file >> rf; file >> n;
			Sphere* sphere = new Sphere(name, x, y, z, i, j, k, r, g, b, a, d, s, rf, n);
			spheres.push_back(sphere);
			std::cout << "Sphere: " << name << " " << x << " " << y << " " << z << " " << i << " " << j << " " << k << " " << r << " " << g << " " << b << " " << a << " " << d << " " << s << " " << rf << " " << n << std::endl;  // Debug print
		}
		else if (command == "LIGHT") {
			// buffer values
			char name[21]; file >> name;
			double x, y, z; file >> x; file >> y; file >> z;
			double r, g, b; file >> r; file >> g; file >> b;
			Light* light = new Light(name, x, y, z, r, g, b);
			lights.push_back(light);
			std::cout << "Light: " << name << " " << x << " " << y << " " << z << " " << r << " " << g << " " << b << std::endl;  // Debug print
		}
		else if (command == "BACK") {
			file >> backgroundColor[0]; file >> backgroundColor[1]; file >> backgroundColor[2];
			std::cout << "Background: " << backgroundColor[0] << " " << backgroundColor[1] << " " << backgroundColor[2] << std::endl;  // Debug print
		}
		else if (command == "AMBIENT") {
			file >> ambientIntensity[0]; file >> ambientIntensity[1]; file >> ambientIntensity[2];
			std::cout << "Ambient: " << ambientIntensity[0] << " " << ambientIntensity[1] << " " << ambientIntensity[2] << std::endl;  // Debug print
		}
		else if (command == "OUTPUT") {
			file >> output_name;
			output_name[sizeof(output_name) - 1] = '\0';
			std::cout << "Output: " << output_name << std::endl;  // Debug print
		}
		else {
			std::cout << "Invalid Input Detected." << std::endl;
		}
	}

	// determines the size of each pixel. Plane is a given space, but resolution determines how that plane is split up
	double scale = resolution[1] / resolution[0];	
	unsigned char *pixels;
	
	pixels = new unsigned char [3 * resolution[0] * resolution[1]]; // 3 values per pixel (r,g,b)

	int k = 0;
	for(int i = 0; i < resolution[1]; i++) {	// for each column of pixels
		for (int j = 0; j < resolution[0]; j++) {	// for each pixel in each row
			double tmin = maxDepth;
			Sphere* nearSphere = nullptr;
			double hitPoint[3];

			Ray* ray = new Ray(j, resolution[1] - 1 - i);
			for (const auto& sphere : spheres) {	// for each sphere
				double t = maxDepth;
				if (sphere->intersect(*ray, t)) {	// if ray intersects with a sphere
					if (t < tmin) {
						tmin = t;
						nearSphere = sphere;	// if it has the closest z value, set it as the sphere to be rendered
						for (int i = 0; i < 3; i++) {
							hitPoint[i] = ray->point[i] * tmin;	// calculate hitpoint/intersection
						}
                	}
				}
			}
			if (nearSphere != nullptr) {
				double pixelColor[3];
            	// sphere was intersected, set pixel color based on closest sphere
				for (int i = 0; i < 3; i++) {
					pixels[k + i] = nearSphere->pixelColor(lights, hitPoint, pixelColor, spheres)[i] * 255;
				}
        	} else {
            // no intersection, set pixel color to background color
				for (int i = 0; i < 3; i++) {
					pixels[k + i] = backgroundColor[i] * 255;
				}
			}
			k = k + 3;
		}
	}

	save_imageP6(resolution[0], resolution[1], output_name, pixels);
	std::cout << "Image rendered.";

	// Clean-up
	for (const auto& sphere : spheres) {
		delete sphere;
	}
	for (const auto& light : lights) {
		delete light;
	}

	file.close();

	return 0;
}