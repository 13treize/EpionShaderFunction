//both function are modified version of https://www.shadertoy.com/view/XtV3z3
float lum(float3 color)
{
    float3 rgb = color;
    return 0.2126 * rgb.r + 0.7152 * rgb.g + 0.0722 * rgb.b;
}