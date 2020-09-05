/*
    Reference materials??
        https://knowledge.autodesk.com/ja/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2015/JPN/Maya/files/Shading-Nodes-htm.html

    Use Support function

    Diffuse

    GammaCorrect

    Luminance
*/
void Diffuse(float3 N, float3 L, float3 C, float3 K, out float3 Out)
{
    float D = dot(N, -L);
    D = max(0, D);
    Out = K * C * D;
}
void Diffuse2(float3 normal, float3 lightvec, float3 diffuse, out float3 Out)
{
    float nDotL = dot(normal, lightvec);
    Out = clamp(nDotL * diffuse, 0.0, 1.0);
}
void GammaCorrect(float3 color,out float3 Out)
{
    Out= pow(color, float3(1.0 / 2.2));
}

void Luminance(float3 color, out float3 Out)
{
    Out = (color.r * 0.3) + (color.g * 0.59) + (color.b * 0.11);
}