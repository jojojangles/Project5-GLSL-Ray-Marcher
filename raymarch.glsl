
// -------------------------
// Distance Functions, from IQ
// -------------------------

float sdPlane( vec3 p )
{
	return p.y;
}

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float sdEllipsoid( in vec3 p, in vec3 r )
{
    return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z);
}

float udRoundBox( vec3 p, vec3 b, float r )
{
  return length(max(abs(p)-b,0.0))-r;
}

float sdTorus( vec3 p, vec2 t )
{
  return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

float sdHexPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
#if 0
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
#else
    float d1 = q.z-h.y;
    float d2 = max((q.x*0.866025+q.y*0.5),q.y)-h.x;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
#endif
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}

float sdTriPrism( vec3 p, vec2 h )
{
    vec3 q = abs(p);
#if 0
    return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
#else
    float d1 = q.z-h.y;
    float d2 = max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
#endif
}

float sdCylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdCone( in vec3 p, in vec3 c )
{
    vec2 q = vec2( length(p.xz), p.y );
    float d1 = -q.y-c.z;
    float d2 = max( dot(q,c.xy), q.y);
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

float sdConeSection( in vec3 p, in float h, in float r1, in float r2 )
{
    float d1 = -p.y - h;
    float q = p.y - h;
    float si = 0.5*(r1-r2)/h;
    float d2 = max( sqrt( dot(p.xz,p.xz)*(1.0-si*si)) + q*si - r2, q );
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

// ------------------------
// Operators
// ------------------------
float opS( float d1, float d2 )
{
    return max(-d2,d1);
}

vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}

vec3 opRep( vec3 p, vec3 c )
{
    return mod(p,c)-0.5*c;
}

vec3 opTwist( vec3 p )
{
    float  c = cos(10.0*p.y+10.0);
    float  s = sin(10.0*p.y+10.0);
    mat2   m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
}

// ------------------------
// Scene Setup
// ------------------------
vec2 map( in vec3 pos )
{
    vec2 res = opU( vec2( opS(sdPlane(pos),sdCylinder(pos, vec2(0.5,0.1))), 1.0 ),
	                vec2( sdSphere(pos-vec3(-1.0,0.15, 0.0), 0.25 ), 50 ) );
    res = opU( res,
               vec2( sdBox(pos-vec3(1.0,.5,-0.5), vec3(0.25,0.1,0.25)), 25) );
    return res;
}

// ------------------------
// Rendery Bits
// ------------------------
const vec3 dirLight = normalize(vec3(-5.0,15.0,-5.0)); //this is TO light, kind of
const float tmin = 1.0;
const float tmax = 20.0;
const float precis = 0.002;
const float eps = 0.0001;
const float k = 128.0; //ao, soft shadows

vec2 castRay( in vec3 ro, in vec3 rd )
{
    float t = tmin;
    float m = -1.0;
    for( int i=0; i<50; i++ )
    {
	    vec2 res = map( ro+rd*t );
        if( res.x<precis || t>tmax ) break;
        t += res.x;
	    m = res.y;
    }

    if( t>tmax ) m=-1.0;
    return vec2( t, m );
}

vec3 getNormal(in vec3 p )
{
    vec2 hxe = map(p + vec3(eps,0,0));
    vec2 hye = map(p + vec3(0,eps,0));
    vec2 hze = map(p + vec3(0,0,eps));

    vec3 n = (1.0/eps) * (map(p).x - vec3(hxe.x,hye.x,hze.x));
    return normalize( -n );
}

float getShadow(in vec3 ro, in vec3 rd) { //ro - point, rd - light direction
    float res = 1.0;
    vec2 ray = castRay(ro,rd);
    float t = ray.x;
    float m = ray.y;
    if(m>-0.5) {
        return 0.0;
    }
    res = min(res,k*m/t);
    return res;
}

float getLighting(in vec3 n, in vec3 p) {
    return dot(n,dirLight) * getShadow(p,dirLight);
}

float getAO(in vec3 p) {
    //separately calculate the 5 'feelers'?
    //totally should have done this in a loop :(
    vec3 n = getNormal(p);
    float f = 0.5 * (eps - castRay(p + n*eps,-n).x);
    f = f + 0.25 * (2.0*eps - castRay(p + n*2.0*eps,-n).x);
    f = f + 0.125 * (3.0*eps - castRay(p + n*3.0*eps,-n).x);
    f = f + 0.0625 * (4.0*eps - castRay(p + n*4.0*eps,-n).x);
    f = f + 0.03125 * (5.0*eps - castRay(p + n*5.0*eps,-n).x);
    return 1.0 - k*f;
}

vec3 render(in vec3 ro, in vec3 rd) {
    //return rd;  // camera ray direction debug view
    vec3 col = vec3(0.1, 0.1, 0.1);
    vec2 res = castRay(ro,rd);
    float t = res.x;
	float m = res.y;
    vec3 p = ro + t*rd;
    if( m>-0.5 ) {
        col = vec3(0.85,0.85,0.85)
            * getLighting(getNormal(p),p);
    }
    return vec3( clamp(col,0.0,1.0) );
}

mat3 setCamera(in vec3 ro, in vec3 ta, float cr) {
    // Starter code from iq's Raymarching Primitives
    // https://www.shadertoy.com/view/Xds3zN

    vec3 cw = normalize(ta - ro);
    vec3 cp = vec3(sin(cr), cos(cr), 0.0);
    vec3 cu = normalize(cross(cw, cp));
    vec3 cv = normalize(cross(cu, cw));
    return mat3(cu, cv, cw);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    // Starter code from iq's Raymarching Primitives
    // https://www.shadertoy.com/view/Xds3zN

    vec2 q = fragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0 * q;
    p.x *= iResolution.x / iResolution.y;
    vec2 mo = iMouse.xy / iResolution.xy;

    float time = 15.0 + iGlobalTime;

    // camera
    vec3 ro = vec3(
            -0.5 + 3.5 * cos(0.1 * time + 6.0 * mo.x),
            1.0 + 2.0 * mo.y,
            0.5 + 3.5 * sin(0.1 * time + 6.0 * mo.x));
    vec3 ta = vec3(-0.5, -0.4, 0.5);

    // camera-to-world transformation
    mat3 ca = setCamera(ro, ta, 0.0);

    // ray direction
    vec3 rd = ca * normalize(vec3(p.xy, 2.0));

    // render
    vec3 col = render(ro, rd);

    col = pow(col, vec3(0.4545));

    fragColor = vec4(col, 1.0);
}
