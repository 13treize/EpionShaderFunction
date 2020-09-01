/*
    Reference materials
        https://thebookofshaders.com
    Hash... https://www.shadertoy.com/view/4djSRW
        exsample hash13 outtype 1(float),intype 3(vec3)
    Use Support function

    Random
        Random1
        Random2
        Random3
    RandomMaze
    Noise
        GradientNoise
        PhasorNoise
        SimpleNoise
        ValueNoise
        Voronoi
            ChebychevVoronoi
            ManhattanVoronoi
            MinkowskiVoronoi
        WaveletNoise
    FBM
        FBM24 -in2,roop4
        FBM28 -in2,roop8
        FBM34 -in3,roop4
        FBM38 -in3,roop8

*/

//Hash
void hash11(float p, out float Out)
{
    p = frac(p * .1031);
    p *= p + 33.33;
    p *= p + p;
    Out = frac(p);
}
void hash12(float2 p, out float Out)
{
    float3 p3 = frac(float3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    Out = frac((p3.x + p3.y) * p3.z);
}
void hash13(float3 p3, out float Out)
{
    p3 = frac(p3 * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    Out = frac((p3.x + p3.y) * p3.z);
}
void hash21(float p, out float2 Out)
{
    float3 p3 = frac(p * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    Out = frac((p3.xx + p3.yz) * p3.zy);
}
void hash22(float2 p, out float2 Out)
{
    float3 p3 = frac(float3(p.xyx) * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    Out = frac((p3.xx + p3.yz) * p3.zy);
}
void hash23(float3 p3, out float2 Out)
{
    p3 = frac(p3 * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    Out = frac((p3.xx + p3.yz) * p3.zy);
}
void hash31(float p, out float3 Out)
{
    float3 p3 = frac(p * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx + 33.33);
    Out = frac((p3.xxy + p3.yzz) * p3.zyx);
}
void hash32(float2 p, out float3 Out)
{
    float3 p3 = frac(float3(p.xyx) * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz + 33.33);
    Out = frac((p3.xxy + p3.yzz) * p3.zyx);
}
void hash33(float3 p3, out float3 Out)
{
    p3 = frac(p3 * float3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz + 33.33);
    Out = frac((p3.xxy + p3.yzz) * p3.zyx);
}
void hash41(float p, out float4 Out)
{
    float4 p4 = frac(p * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    Out = frac((p4.xxyz + p4.yzzw) * p4.zywx);
}
void hash42(float2 p, out float4 Out)
{
    float4 p4 = frac(float4(p.xyxy) * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    Out = frac((p4.xxyz + p4.yzzw) * p4.zywx);
}
void hash43(float3 p, out float4 Out)
{
    float4 p4 = frac(float4(p.xyzx) * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    Out = frac((p4.xxyz + p4.yzzw) * p4.zywx);
}
void hash44(float4 p4, out float4 Out)
{
    p4 = frac(p4 * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    Out = frac((p4.xxyz + p4.yzzw) * p4.zywx);
}

/*
    use support function
    Naming to avoid name calling
*/
float NoiseMaxAbs(float2 p)
{
    return max(abs(p.x), abs(p.y));
}
float NoiseSumAbs(float2 p)
{
    return abs(p.x) + abs(p.y);
}
float NoiseLpNorm(float2 p, float n)
{
    float2 t = pow(abs(p), n);
    return pow(t.x + t.y, 1.0 / n);
}
float2x2 NoiseRot(float a)
{
    float s = sin(a);
    float c = cos(a);
    return float2x2(c, -s, s, c);
}

//Random
float random_seed1(float val)
{
    return frac(sin(val) *  43758.5453123);
}
float random_seed2(float2 uv)
{
    return frac(sin(dot(uv.xy, float2(12.9898, 78.233))) * 43758.5453123);
}

float random_seed3(float3 p)
{
    return frac(sin(dot(p.xyz, float3(12.9898, 78.233,54.641))) * 43758.5453123);
}
void Random(float2 uv, float2 scale, out float iOut, out float fOut)
{
    float2 p = uv * scale;
    iOut = random_seed2(floor(p));
    fOut = random_seed2(frac(p));
}

//RandomMaze
float2 truchetPattern(float2 uv, float index)
{
    index = frac(((index - 0.5) * 2.0));
    if (index > 0.75)
    {
        uv = 1.0 - uv;
    }
    else if (index > 0.5)
    {
        uv = float2(1.0 - uv.x, uv.y);
    }
    else if (index > 0.25)
    {
        uv = 1.0 - float2(1.0 - uv.x, uv.y);
    }
    return uv;
}

void RandomMaze(float2 uv, float2 scale, out float Line, out float Circle, out float Truchet)
{
    float I, F;
    Random(uv, scale, I, F);
    float2 p = uv * scale;
    float2 tile = truchetPattern(frac(p), I);
    Line = smoothstep(tile.x - 0.3, tile.x, tile.y) - smoothstep(tile.x, tile.x + 0.3, tile.y);
    Circle = (step(length(tile), 0.6) - step(length(tile), 0.4)) + (step(length(tile - 1.0), 0.6) - step(length(tile - 1.0), 0.4));
    Truchet = step(tile.x, tile.y);
}


//GradientNoise
float2 gradientNoise_dir(float2 p)
{
    p = p % 289;
    float x = (34 * p.x + 1) * p.x % 289 + p.y;
    x = (34 * x + 1) * x % 289;
    x = frac(x / 41) * 2 - 1;
    return normalize(float2(x - floor(x + 0.5), abs(x) - 0.5));
}
float gradient_noise(float2 p)
{
    float2 ip = floor(p);
    float2 fp = frac(p);
    float d00 = dot(gradientNoise_dir(ip), fp);
    float d01 = dot(gradientNoise_dir(ip + float2(0, 1)), fp - float2(0, 1));
    float d10 = dot(gradientNoise_dir(ip + float2(1, 0)), fp - float2(1, 0));
    float d11 = dot(gradientNoise_dir(ip + float2(1, 1)), fp - float2(1, 1));
    fp = fp * fp * fp * (fp * (fp * 6 - 15) + 10);
    return lerp(lerp(d00, d01, fp.y), lerp(d10, d11, fp.y), fp.x);
}
void GradientNoise(float2 UV, float Scale, out float Out)
{
    Out = gradient_noise(UV * Scale) + 0.5;
}
//PhasorNoise    èàóùÇ™èdÇ¢
int morton(int x, int y)
{
    int z = 0;
    for (int i = 0; i < 32 * 4; i++)
    {
        z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
    }
    return z;
}
float2 phasor(float2 x, float f, float b, float o, float phi)
{
    float a = exp(-3.14 * (b * b) * ((x.x * x.x) + (x.y * x.y))),
          s = sin(2.0 * 3.14 * f * (x.x * cos(o) + x.y * sin(o)) + phi),
          c = cos(2.0 * 3.14 * f * (x.x * cos(o) + x.y * sin(o)) + phi);
    return float2(a * c, a * s);
}
float uni_0_1(inout int xx)
{
    int N = 15487469;
    xx *= 3039177861;
    xx %= N;
    return float(xx) / float(N);
}
float uni(inout int xx, float Min, float Max)
{
    float uni_0_1_ = uni_0_1(xx);
    return Min + (uni_0_1_ * (Max - Min));
}
float2 cell(int2 ij, float2 uv, float _kr, float f, float b, float o)
{
    int s = morton(ij.x, ij.y) + 333;
    s = s == 0 ? 1 : s + 2;
    int x = s, impulse = 0, nImpulse = 16;
    float cellsz = 2.0 * _kr;
    float2 noise;
    while (impulse <= nImpulse)
    {
        float2 impulse_centre = float2(uni_0_1(x), uni_0_1(x));
        float2 d = (uv - impulse_centre) * cellsz;
        float rp = uni(x, 0.0, 2.0 * 3.14);
        noise += phasor(d, f, b, o, rp);
        impulse++;
    }
    return noise;
}
float2 eval_noise(float2 UV, float kr, float f, float b, float o)
{

    float cellsz = kr * 2.0;
    float2 ij_ = UV / cellsz;
    int2 ij = int2(int(ij_.x), int(ij_.y));
    float2 fij = ij_ - float2(float(ij.x), float(ij.y));
    float2 noise;
    for (int j = -2; j <= 2; j++)
    {
        for (int i = -2; i <= 2; i++)
        {
            int2 nij = int2(i, j);
            noise += cell(ij + nij, fij - float2(float(nij.x), float(nij.y)), kr, f, b, o);
        }
    }
    return noise;
}
void PhasorNoise(float2 uv, float f_, float b_, float o_, out float a, out float b, out float c, out float d)
{
    float2 dir = float2(cos(o_), sin(o_));
    float kr = sqrt(-log(0.05) / 3.14) / b_;
    float2 phasorNoise = eval_noise(uv, kr, f_, b_, o_) * 0.3;
    float phi = atan2(phasorNoise.y, phasorNoise.x);
    a = phi;
    b = sin(phi) + 0.5;
    c = length(phasorNoise);
    d = frac(phi / (2.0 * 3.14) - f_ * dot(uv, dir));
}
//  SimpleNoise
inline float noise_randomValue(float2 uv)
{
    return frac(sin(dot(uv, float2(12.9898, 78.233))) * 43758.5453);
}
inline float noise_interpolate(float a, float b, float t)
{
    return (1.0 - t) * a + (t * b);
}
inline float valueNoise(float2 uv)
{
    float2 i = floor(uv);
    float2 f = frac(uv);
    f = f * f * (3.0 - 2.0 * f);
    uv = abs(frac(uv) - 0.5);
    float2 c0 = i + float2(0.0, 0.0);
    float2 c1 = i + float2(1.0, 0.0);
    float2 c2 = i + float2(0.0, 1.0);
    float2 c3 = i + float2(1.0, 1.0);
    float r0 = noise_randomValue(c0);
    float r1 = noise_randomValue(c1);
    float r2 = noise_randomValue(c2);
    float r3 = noise_randomValue(c3);
    float bottomOfGrid = noise_interpolate(r0, r1, f.x);
    float topOfGrid = noise_interpolate(r2, r3, f.x);
    float t = noise_interpolate(bottomOfGrid, topOfGrid, f.y);
    return t;
}
void SimpleNoise(float2 UV, float Scale, out float Out)
{
    float t = 0.0;
    float freq = pow(2.0, float(0));
    float amp = pow(0.5, float(3 - 0));
    t += valueNoise(float2(UV.x * Scale / freq, UV.y * Scale / freq)) * amp;
    freq = pow(2.0, float(1));
    amp = pow(0.5, float(3 - 1));
    t += valueNoise(float2(UV.x * Scale / freq, UV.y * Scale / freq)) * amp;
    freq = pow(2.0, float(2));
    amp = pow(0.5, float(3 - 2));
    t += valueNoise(float2(UV.x * Scale / freq, UV.y * Scale / freq)) * amp;
    Out = t;
};
//  Voronoi
float2 voronoi_noise_randomVector(float2 UV, float2 offset)
{
    float2x2 m = float2x2(15.27, 47.63, 99.41, 89.98);
    UV = frac(sin(mul(UV, m)) * 46839.32);
    return float2(sin(UV.y * +offset.x) * 0.5 + 0.5, cos(UV.x * offset.y) * 0.5 + 0.5);
}
//  Voronoi
void Voronoi(float2 UV, float2 AngleOffset, float2 CellDensity, out float Out, out float Cells, out float Lines, out float Points)
{
    float2 g = floor(UV * CellDensity);
    float2 f = frac(UV * CellDensity);
    float res = 8.0;
    float md = 8.0;
    float2 mr;
    for (int y = -1; y <= 1; y++)
    {
        for (int x = -1; x <= 1; x++)
        {
            float2 lattice = float2(x, y);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = dot(r, r);
            if (d < res)
            {
                res = d;
                mr = r;
            }
        }
    }
    res = 8.0;
    for (int yy = -1; yy <= 1; yy++)
    {
        for (int xx = -1; xx <= 1; xx++)
        {
            float2 lattice = float2(xx, yy);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = dot(r, r);
            if (d < res)
            {
                res = d;
                Out = res;
                Cells = offset.x;
            }
            if (dot(mr - r, mr - r) > 0.00001)
            {
                md = min(md, dot(0.5 * (mr + r), normalize(r - mr)));
            }
        }
    }
    Lines = lerp(1.0, 0.0, smoothstep(0.03, 0.06, md));
    Points = 1.0 - smoothstep(0.0, 0.1, res);
}
//  ChebychevVoronoi
void ChebychevVoronoi(float2 UV, float2 AngleOffset, float2 CellDensity, out float Out, out float Cells, out float Lines, out float Points)
{
    float2 g = floor(UV * CellDensity);
    float2 f = frac(UV * CellDensity);
    float res = 8.0;
    float md = 8.0;
    float2 mr;
    for (int y = -1; y <= 1; y++)
    {
        for (int x = -1; x <= 1; x++)
        {
            float2 lattice = float2(x, y);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = NoiseMaxAbs(r);
            if (d < res)
            {
                res = d;
                mr = r;
            }
        }
    }
    res = 8.0;
    for (int yy = -1; yy <= 1; yy++)
    {
        for (int xx = -1; xx <= 1; xx++)
        {
            float2 lattice = float2(xx, yy);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = NoiseMaxAbs(r);
            if (d < res)
            {
                res = d;
                Out = res;
                Cells = offset.x;
            }
            if (dot(mr - r, mr - r) > 0.00001)
            {
                md = min(md, dot(0.5 * (mr + r), normalize(r - mr)));
            }
        }
    }
    Lines = lerp(1.0, 0.0, smoothstep(0.03, 0.06, md));
    Points = 1.0 - smoothstep(0.0, 0.1, res);
}
//  ManhattanVoronoi
void ManhattanVoronoi(float2 UV, float2 AngleOffset, float2 CellDensity, out float Out, out float Cells, out float Lines, out float Points)
{
    float2 g = floor(UV * CellDensity);
    float2 f = frac(UV * CellDensity);
    float res = 8.0;
    float md = 8.0;
    float2 mr;
    for (int y = -1; y <= 1; y++)
    {
        for (int x = -1; x <= 1; x++)
        {
            float2 lattice = float2(x, y);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = NoiseSumAbs(r);
            if (d < res)
            {
                res = d;
                mr = r;
            }
        }
    }
    res = 8.0;
    for (int yy = -1; yy <= 1; yy++)
    {
        for (int xx = -1; xx <= 1; xx++)
        {
            float2 lattice = float2(xx, yy);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = NoiseSumAbs(r);
            if (d < res)
            {
                res = d;
                Out = res;
                Cells = offset.x;
            }
            if (dot(mr - r, mr - r) > 0.00001)
            {
                md = min(md, dot(0.5 * (mr + r), normalize(r - mr)));
            }
        }
    }
    Lines = lerp(1.0, 0.0, smoothstep(0.03, 0.06, md));
    Points = 1.0 - smoothstep(0.0, 0.1, res);
}
// MinkowskiVoronoi
void MinkowskiVoronoi(float2 UV, float2 AngleOffset, float2 CellDensity, float p, out float Out, out float Cells, out float Lines, out float Points)
{
    float2 g = floor(UV * CellDensity);
    float2 f = frac(UV * CellDensity);
    float res = 8.0;
    float md = 8.0;
    float2 mr;
    for (int y = -1; y <= 1; y++)
    {
        for (int x = -1; x <= 1; x++)
        {
            float2 lattice = float2(x, y);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = NoiseLpNorm(r, p);
            if (d < res)
            {
                res = d;
                mr = r;
            }
        }
    }
    res = 8.0;
    for (int yy = -1; yy <= 1; yy++)
    {
        for (int xx = -1; xx <= 1; xx++)
        {
            float2 lattice = float2(xx, yy);
            float2 offset = voronoi_noise_randomVector(lattice + g, AngleOffset);
            float2 r = lattice + offset - f;
            float d = NoiseLpNorm(r, p);
            if (d < res)
            {
                res = d;
                Out = res;
                Cells = offset.x;
            }
            if (dot(mr - r, mr - r) > 0.00001)
            {
                md = min(md, dot(0.5 * (mr + r), normalize(r - mr)));
            }
        }
    }
    Lines = lerp(1.0, 0.0, smoothstep(0.03, 0.06, md));
    Points = 1.0 - smoothstep(0.0, 0.1, res);
}

//WaveletNoise
float GetWavelet(float2 p, float angle, float z, float scale)
{
    p *= scale;
    float2 id = floor(p);
    float n;
    hash12(id, n);
    p = frac(p) - 0.5;
    p = mul(p, NoiseRot(n * angle));
    float d = sin(p.x * 10.0 + z);
    d *= smoothstep(0.25, 0.0, dot(p, p));
    float Out = d / scale;
    return Out;
}

void WaveletNoise(float2 p, float angle, float phase, float seed_scale, float scale_factor, out float Out)
{
    float d = 0.0;
    float scale = seed_scale;
    float mag = 0.;
    for (float i = 0.; i < 4.; i++)
    {
        d += GetWavelet(p, angle, phase, scale);
        p = mul(p, float2x2(.54, -.84, .84, .54)) + i;
        mag += 1. / scale;
        scale *= scale_factor;
    }
    d /= mag;
    Out = d * .5 + .5;
}
//FBM
float fbm_noise2(float2 p)
{
    float2 ip = floor(p);
    float2 u = frac(p);
    u = u * u * (3.0 - 2.0 * u);
    float res = lerp(lerp(random_seed2(ip), random_seed2(ip + float2(1.0, 0.0)), u.x),lerp(random_seed2(ip + float2(0.0, 1.0)), random_seed2(ip + float2(1.0, 1.0)), u.x), u.y);
    return res * res;
}

void FBM24(float2 uv, float amplitude, float frequency, out float Out)
{
    float2 p = (uv * 2.0 - 1.0);
    float result = 0.;
    float amplitude2 = amplitude;
    float frequency2 = frequency;
    for (int i = 0; i < 4; i++)
    {
        result += fbm_noise2(p * frequency2) * amplitude2;
        amplitude2 *= .5;
        frequency2 *= 2.;
    }
    Out = result;
}
void FBM28(float2 uv, float amplitude, float frequency, out float Out)
{
    float2 p = (uv * 2.0 - 1.0);
    float result = 0.;
    float amplitude2 = amplitude;
    float frequency2 = frequency;
    for (int i = 0; i < 8; i++)
    {
        result += fbm_noise2(p * frequency2) * amplitude2;
        amplitude2 *= .5;
        frequency2 *= 2.;
    }
    Out = result;
}


//void FBM38(float3 p, float amplitude, float frequency, out float Out)
//{
//    float3 pp = (p * 2.0 - 1.0);
//    float result = 0.;
//    float amplitude2 = amplitude;
//    float frequency2 = frequency;
//    for (int i = 0; i < 8; i++)
//    {
//        result += fbm_noise(p * frequency2) * amplitude2;
//        amplitude2 *= .5;
//        frequency2 *= 2.;
//    }
//    Out = result;
//}