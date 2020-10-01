/*
    Reference materials

    Use Support function


*/
float bias(float a, float b)
{
    return pow(a, log(b) / log(0.5));
}


void Checkerboard(float2 UV, float3 ColorA, float3 ColorB, float2 Frequency, out float3 Out)
{
    UV = (UV.xy + 0.5) * Frequency;
    float4 derivatives = float4(ddx(UV), ddy(UV));
    float2 duv_length = sqrt(float2(dot(derivatives.xz, derivatives.xz), dot(derivatives.yw, derivatives.yw)));
    float width = 1.0;
    float2 distance3 = 4.0 * abs(frac(UV + 0.25) - 0.5) - width;
    float2 scale = 0.35 / duv_length.xy;
    float freqLimiter = sqrt(clamp(1.1f - max(duv_length.x, duv_length.y), 0.0, 1.0));
    float2 vector_alpha = clamp(distance3 * scale.xy, -1.0, 1.0);
    float alpha = saturate(0.5f + 0.5f * vector_alpha.x * vector_alpha.y * freqLimiter);
    Out = lerp(ColorA, ColorB, alpha.xxx);
}

void Ellipse(float2 UV, float Width, float Height, out float Out)
{
    float d = length((UV * 2 - 1) / float2(Width, Height));
    Out = saturate((1 - d) / fwidth(d));
}

void Gear(float2 UV, float tooths, float radius, float tlength, float tshape, float intooths, float inradius, out float Out)
{
    float2 p = UV * 2.0 - 1.0;
    float col = 1.0;
    float gear = max(radius, 0.15) + tlength * clamp((cos(atan2(p.x, p.y) * tooths)), 0.0, clamp(tshape, 0.08, 1.0));
    gear = min(min(gear, step(0.03, length(p))), step(clamp(inradius, 0.1, radius) * (sign(cos((atan2(p.x, p.y) * intooths)))), length(p)));
    col *= min(step(gear, length(p)), step(0.05, length(p)));
    Out = 1.0 - col;
}
void Hexagon(float2 UV, float Scale, out float Out, out float2 Pos, out float2 oUV, out float2 Index)
{
    float2 p = UV * Scale;
    p.x *= 1.15470053838;
    float isTwo = frac(floor(p.x) / 2.0) * 2.0;
    float isOne = 1.0 - isTwo;
    p.y += isTwo * 0.5;
    float2 rectUV = frac(p);
    float2 grid = floor(p);
    p = frac(p) - 0.5;
    float2 s = sign(p);
    p = abs(p);
    Out = abs(max(p.x * 1.5 + p.y, p.y * 2.0) - 1.0);
    float isInHex = step(p.x * 1.5 + p.y, 1.0);
    float isOutHex = 1.0 - isInHex;
    float2 grid2 = float2(0, 0);
    grid2 = lerp(float2(s.x, +step(0.0, s.y)), float2(s.x, -step(s.y, 0.0)), isTwo) * isOutHex;
    Index = grid + grid2;
    Pos = Index / Scale;
    oUV = lerp(rectUV, rectUV - s * float2(1.0, 0.5), isOutHex);
}

void Polygon(float2 UV, float Sides, float Width, float Height, out float Out)
{
    float pi = 3.14159265359;
    float aWidth = Width * cos(pi / Sides);
    float aHeight = Height * cos(pi / Sides);
    float2 uv = (UV * 2 - 1) / float2(aWidth, aHeight);
    uv.y *= -1;
    float pCoord = atan2(uv.x, uv.y);
    float r = 2 * pi / Sides;
    float distance = cos(floor(0.5 + pCoord / r) * r - pCoord) * length(uv);
    Out = saturate((1 - distance) / fwidth(distance));
};

void Ripple(float2 UV, float Width, float Height, float2 Center, float Scale, out float Out)
{
    float2 p = (UV * 2.0 - 1.0) / float2(Width, Height);
    Out = sin(Scale * distance(p, Center));
}

void RoundedRectangle(float2 UV, float Width, float Height, float Radius, out float Out)
{
    Radius = max(min(min(abs(Radius * 2), abs(Width)), abs(Height)), 1e-5);
    float2 uv = abs(UV * 2 - 1) - float2(Width, Height) + Radius;
    float d = length(max(0, uv)) / Radius;
    Out = saturate((1 - d) / fwidth(d));
}