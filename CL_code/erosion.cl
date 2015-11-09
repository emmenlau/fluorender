#define KX 3
#define KY 3
#define KZ 3
const sampler_t samp =
	CLK_NORMALIZED_COORDS_FALSE|
	CLK_ADDRESS_CLAMP_TO_EDGE|
	CLK_FILTER_LINEAR;
__kernel void kernel_main(
	read_only image3d_t data,
	__global unsigned char* result,
	unsigned int x,
	unsigned int y,
	unsigned int z)
{
	int4 coord = (int4)(get_global_id(0),
		get_global_id(1), get_global_id(2), 1);
	int4 kc;
	float4 dvalue;
	float rvalue = 1.0;
	float r = 0.5;
	int i, j, k;
	for (i=0; i<KX; ++i)
	for (j=0; j<KY; ++j)
	for (k=0; k<KZ; ++k)
	{
		float l = sqrt((float)((i-KX/2)*(i-KX/2)+
				(j-KY/2)*(j-KY/2)+
				(k-KZ/2)*(k-KZ)/2));
		if (l == 0.0)
			continue;
		kc = (int4)(coord.x+(i-KX/2)*r/l,
				coord.y+(j-KY/2)*r/l,
				coord.z+(k-KZ/2)*r/l, 1);
		dvalue = read_imagef(data, samp, kc);
		rvalue = min(rvalue, dvalue.x);
	}
	unsigned int index = x*y*coord.z + x*coord.y + coord.x;
	result[index] = rvalue*255.0;
}