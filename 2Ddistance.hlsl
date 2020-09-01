/*
    Use Support function
    Box2D
    Circle2D
    RoundedBox2D
*/
void Box2D(float2 uv, float2 b, out float Out)
{
    float2 d = abs(uv) - b;
    Out = length(max(d, 0.0)) + min(max(d.x, d.y), 0.0);
}

void Circle2D(float2 uv, float r,out float Out)
{
    Out= length(uv) - r;
}

void RoundedBox2D(float2 uv, float2 b, float4 r, out float Out)
{
    r.xy = (uv.x > 0.0) ? r.xy : r.zw;
    r.x = (uv.y > 0.0) ? r.x : r.y;
    float2 q = abs(uv) - b + r.x;
    Out= min(max(q.x, q.y), 0.0) + length(max(q, 0.0)) - r.x;
}
