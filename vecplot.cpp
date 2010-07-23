// Compile with g++ -O2 *.cpp kmpp/*.cpp -o vecplot

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>

typedef struct{
	int n_arrows;
	int kmeans_attempts;
	int use_all_vectors;
} options;

typedef struct{
	int n;
	double *xy, *f;
	double max_len;
	double bbox[4]; // left, right, bottom, top
} data;

static double pythag2(double x, double y){
	double ax = fabs(x);
	double ay = fabs(y);
	if(0 == ax && 0 == ay){ return 0; }
	if(ax > ay){
		double r = ay/ax;
		return ax*sqrt(1+r*r);
	}else{
		double r = ax/ay;
		return ay*sqrt(1+r*r);
	}
}

void data_init(data *d){
	if(NULL == d){ return; }
	d->n = 0;
	d->xy = NULL;
	d->f = NULL;
	d->max_len = 0;
	d->bbox[0] = DBL_MAX;
	d->bbox[1] = -DBL_MAX;
	d->bbox[2] = DBL_MAX;
	d->bbox[3] = -DBL_MAX;
}
void data_destroy(data *d){
	if(NULL == d){ return; }
	if(NULL != d->xy){ free(d->xy); }
}

int data_read(data *d, const char *filename){
	FILE *fp;
	int i;
	char line[1024];
	int ncap = 256;
	int line_count = 0, count;
	double x, y, fx, fy, fa;
	double *temp;
	
	if(NULL == d){ return -1; }
	if(NULL == filename){ return -2; }
	
	d->n = 0;
	if(NULL != d->xy){ free(d->xy); }
	temp = (double*)malloc(sizeof(double)*4*ncap);
	
	fp = fopen(filename, "rt");
	if(NULL == fp){ return -3; }
	
	while(fgets(line, sizeof(line), fp) != NULL){ ++line_count;
		if('#' == line[0]){ continue; }
		
		count = sscanf(line, "%lf %lf %lf %lf", &x, &y, &fx, &fy);
		if(0 < count && count < 4){
			fprintf(stderr, "Expected 4 values on line %d but only got %d\n", line_count, count);
			goto error;
		}
		if(count < 1){ continue; }

		if(d->n >= ncap){
			ncap *= 2;
			temp = (double*)realloc(temp, sizeof(double)*4*ncap);
		}
		
		temp[4*(d->n)+0] = x;
		temp[4*(d->n)+1] = y;
		temp[4*(d->n)+2] = fx;
		temp[4*(d->n)+3] = fy;
		d->n++;
		if(x < d->bbox[0]){ d->bbox[0] = x; }
		if(x > d->bbox[1]){ d->bbox[1] = x; }
		if(y < d->bbox[2]){ d->bbox[2] = y; }
		if(y > d->bbox[3]){ d->bbox[3] = y; }
		fa = pythag2(fx, fy);
		if(fa > d->max_len){ d->max_len = fa; }
	}
	// Enlarge bounding box slightly
	d->bbox[0] -= DBL_EPSILON*d->bbox[0];
	d->bbox[1] += DBL_EPSILON*d->bbox[1];
	d->bbox[2] -= DBL_EPSILON*d->bbox[2];
	d->bbox[3] += DBL_EPSILON*d->bbox[3];
	
	// Copy data into final arrays
	d->xy = (double*)malloc(sizeof(double)*4*d->n);
	d->f = d->xy + 2*d->n;
	for(i = 0; i < d->n; ++i){
		d->xy[2*i+0] = temp[4*i+0];
		d->xy[2*i+1] = temp[4*i+1];
		d->f[2*i+0] = temp[4*i+2];
		d->f[2*i+1] = temp[4*i+3];
	}
	free(temp);
error:
	fclose(fp);
}

int gen_plot_kmeans(data *raw, int n, data *plot, int attempts){
	int *which;
	int i,j,k;
	double *work;
	
	extern double RunKMeansPlusPlus(int n, int k, int d, double *points, int attempts,
                 double *centers, int *assignments);

	plot->n = n;
	plot->xy = (double*)malloc(sizeof(double)*4*n);
	plot->f = plot->xy + 2*n;
	which = (int*)malloc(sizeof(int)*raw->n);
	work = (double*)malloc(sizeof(double)*2*n);
	for(i = 0; i < 2*n; ++i){
		work[i] = 0;
		plot->f[i] = 0;
	}
	RunKMeansPlusPlus(raw->n, n, 2, raw->xy, attempts, plot->xy, which);
	for(i = 0; i < raw->n; ++i){
		int c = which[i];
		double d = pythag2(plot->xy[2*c+0] - raw->xy[2*i+0], plot->xy[2*c+1] - raw->xy[2*i+1]);
		if(d > work[c]){ work[c] = d; }
	}
	
	for(i = 0; i < raw->n; ++i){
		int c = which[i];
		double d = pythag2(plot->xy[2*c+0] - raw->xy[2*i+0], plot->xy[2*c+1] - raw->xy[2*i+1]);
		if(0 == work[c]){
			d = 1.0;
		}else{
			d = 1.0 - d/work[c];
		}
		plot->f[2*c+0] += d*raw->f[2*i+0];
		plot->f[2*c+1] += d*raw->f[2*i+1];
		work[n+c] += d;
	}
	
	plot->max_len = 0;
	for(i = 0; i < n; ++i){
		if(work[n+i] > 0){
			plot->f[2*i+0] /= work[n+i];
			plot->f[2*i+1] /= work[n+i];
			double d = pythag2(plot->f[2*i+0], plot->f[2*i+1]);
			if(d > plot->max_len){ plot->max_len = d; }
			if(plot->xy[2*i+0] < plot->bbox[0]){ plot->bbox[0] = plot->xy[2*i+0]; }
			if(plot->xy[2*i+0] > plot->bbox[1]){ plot->bbox[1] = plot->xy[2*i+0]; }
			if(plot->xy[2*i+1] < plot->bbox[2]){ plot->bbox[2] = plot->xy[2*i+1]; }
			if(plot->xy[2*i+1] > plot->bbox[3]){ plot->bbox[3] = plot->xy[2*i+1]; }
		}
		//fprintf(stderr, "%f %f %f %f\n", plot->xy[2*i+0], plot->xy[2*i+1], plot->f[2*i+0], plot->f[2*i+1]);
	}
	free(work);
	free(which);
	
	return 0;
}

void map_bb_aff(const double from[4], const double to[4], double p[2]){
	p[0] = to[0] + (p[0]-from[0])/(from[1]-from[0]) * (to[1]-to[0]);
	p[1] = to[2] + (p[1]-from[2])/(from[3]-from[2]) * (to[3]-to[2]);
}
void map_bb_lin(const double from[4], const double to[4], double p[2]){
	p[0] = p[0]/(from[1]-from[0]) * (to[1]-to[0]);
	p[1] = p[1]/(from[3]-from[2]) * (to[3]-to[2]);
}

double round_to_nice(double x, int *mant, int *base, int *power){
	// if x = a*10^p, where a is [1, 10)
	// then log10(x) = p + log10(a)
	double lx = log10(x);
	double a = fmod(lx, 1.0);
	int p = (int)lx;
	// round a to nearest 1, 2, or 5
	double r;
	if(a < 0.15051499783199057){ // halfway between 1 and 2
		r = 1;
	}else if(a < 0.5){ // halfway between 2 and 5
		r = 2;
	}else if(a < 0.8494850021680094){ // halfwat between 5 and 10
		r = 5;
	}else{
		r = 10;
	}
	*mant = r;
	*base = 10;
	*power = p;
	return r*pow(10.,p);
}
/*
void get_ticks(double a, double b, double *tick, int nticks){
	int m,b,p;
	if(a <= 0 && b >= 0){ // includes zero
		double m = (-a > b) ? -a : b;
		double t = round_to_nice(0.33*m, &m,&b,&p);
		
	}else{
		double d = b-a;
		double avg = 0.5*(b+a);
		if(d < 0.01*avg){
		}else{
			double t = round_to_nice(0.25*m, &m,&b,&p);
			if(a >= 0){ // positive
				d = b-a;
			}else{ // negative
			}
		}
	}
}
*/
int output_plot(data *plot){
	FILE *f = stdout;
	int i, j;
	double norm;
	double p[2], q[2];
	
	const double bbox[4] = {
		-3,3, // x range
		-3,3  // y range
	};
	double scale[2] = {72, 72};
	fprintf(f, "%f %f scale\n", scale[0], scale[1]);
	fprintf(f, "%f setlinewidth\n", 2./scale[0]);
	fprintf(f, "%f %f translate\n", 8.5*0.5*72/scale[0], 11*0.5*72/scale[1]);
	
	fprintf(f,
		"/arrow{\n"
		"gsave\n"
		"5 3 roll translate\n"
		"3 1 roll exch atan rotate\n"
		"dup scale\n"
		"\n"
		"newpath\n"
		"1.00  0.00 moveto\n"
		"0.62  0.19 lineto\n"
		"0.62  0.07 lineto\n"
		"0.00  0.07 lineto\n"
		"0.00 -0.07 lineto\n"
		"0.62 -0.07 lineto\n"
		"0.62 -0.19 lineto\n"
		"closepath stroke\n"
		"\n"
		"grestore\n"
		"} bind def\n");
	
	// find closest pair of points
	double min_spacing = DBL_MAX;
	for(i = 0; i < plot->n; ++i){
		for(j = i+1; j < plot->n; ++j){
			double d = pythag2(plot->xy[2*i+0]-plot->xy[2*j+0], plot->xy[2*i+1]-plot->xy[2*j+1]);
			if(d < min_spacing){ min_spacing = d; }
		}
	}
	
	norm = min_spacing / plot->max_len;

	for(i = 0; i < plot->n; ++i){
		double vec[2] = {
			plot->f[2*i+0] * norm,
			plot->f[2*i+1] * norm };
		double base[2] = {
			plot->xy[2*i+0] - 0.5*vec[0],
			plot->xy[2*i+1] - 0.5*vec[1] };
		map_bb_aff(plot->bbox, bbox, base);
		map_bb_lin(plot->bbox, bbox, vec);
		double len = pythag2(vec[0], vec[1]);
		if(len < 100*DBL_EPSILON){ continue; }
		fprintf(f, "%f %f %f %f %f arrow\n", base[0], base[1], vec[0], vec[1], len);
	}
	
	// Draw frame and axes
	p[0] = plot->bbox[0];
	p[1] = plot->bbox[2];
	q[0] = plot->bbox[1];
	q[1] = plot->bbox[3];
	map_bb_aff(plot->bbox, bbox, p);
	map_bb_aff(plot->bbox, bbox, q);
	fprintf(f, "newpath %f %f moveto %f %f lineto %f %f lineto %f %f lineto closepath stroke\n",
		p[0], p[1],
		q[0], p[1],
		q[0], q[1],
		p[0], q[1]
		);
	/*
	// Draw ticks
	for(i = 0; i < 2; ++i){ // which dimension
		double ticks[8];
		get_ticks(plot->bbox[2*i+0], plot->bbox[2*i+1], ticks, 8);
	}
	*/
	fprintf(f, "showpage\n");
	return 0;
}

int data_cull_zeros(data *d, int n_arrows){
	int i, c = 0;
	int ret = 0;
	int *flag = (int*)malloc(sizeof(int)*d->n);
	for(i = 0; i < d->n; ++i){
		flag[i] = 0;
		if(pythag2(d->f[2*i+0], d->f[2*i+1]) < DBL_EPSILON*d->max_len){
			flag[i] = 1;
			++c;
		}
	}
	if(d->n - c <= n_arrows){
		// Remove all the zeros
		c = 0;
		for(i = 0; i < d->n; ++i){
			if(!flag[i] && c != i){
				d->xy[2*c+0] = d->xy[2*i+0];
				d->xy[2*c+1] = d->xy[2*i+1];
				d->f[2*c+0] = d->f[2*i+0];
				d->f[2*c+1] = d->f[2*i+1];
				++c;
			}
		}
		d->n = c;
		ret = 1;
	}
	free(flag);
	return ret;
}

void usage(){
	fprintf(stderr, "Usage: vecplot [-a] [-n num-arrows] [-k k-means attempts] file file ...\n");
}

int main(int argc, char *argv[]){
	options opts;
	opts.n_arrows = 100;
	opts.kmeans_attempts = 1;
	opts.use_all_vectors = 0;
	
	int c, i;
	opterr = 0;
	while((c = getopt(argc, argv, "ak:n:")) != -1){
		switch(c){
		case 'a':
			opts.use_all_vectors = 1;
			break;
		case 'k':
			opts.kmeans_attempts = atoi(optarg);
			if(opts.kmeans_attempts < 1){
				opts.kmeans_attempts = 1;
			}
			break;
		case 'n':
			opts.n_arrows = atoi(optarg);
			if(opts.n_arrows < 1){
				opts.n_arrows = 1;
			}
			break;
		case '?':
			if('k' == optopt || 'n' == optopt){
				fprintf(stderr, "Option -%c requires an argument.\n", optopt);
			}else if(isprint(optopt)){
				fprintf(stderr, "Unknown option `-%c'.\n", optopt);
			}else{
				fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
			}
		default:
			usage();
			return EXIT_FAILURE;
		}
	}
	
	if(optind >= argc){ usage(); return EXIT_FAILURE; }
	for(i = optind; i < argc; i++){
		data d;
		data_init(&d);
		data_read(&d, argv[i]);
		
		if(opts.use_all_vectors || d.n <= opts.n_arrows || data_cull_zeros(&d, opts.n_arrows)){
			output_plot(&d);
		}else{
			data plot_data;
			data_init(&plot_data);
			gen_plot_kmeans(&d, opts.n_arrows, &plot_data, opts.kmeans_attempts);
			output_plot(&plot_data);
			data_destroy(&plot_data);
		}
		
		data_destroy(&d);
	}
	
	return EXIT_SUCCESS;
}
