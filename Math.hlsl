/*
    Use Support function
    
*/
float2x2 Rotate2d(float angle)
{
    return float2x2(cos(angle), sin(angle), -sin(angle), cos(angle));
}