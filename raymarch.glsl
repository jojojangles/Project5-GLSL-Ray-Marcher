
// -------------------------
// Distance Functions, from IQ
// - returns the distance to the object
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


float hash( float n ) //util for noise
{
    return fract(sin(n)*43758.5453);
}

float lerp(float a, float b, float t) //util for noise
{
    return (1.0-t)*a + t*b;
}

const float noiseScale = 5.0;

float noise( vec3 p ) //make noise a distance function?
{
    // range -1.0f -> 1.0f
    vec3 pf = floor(p*noiseScale);
    vec3 f = fract(p*noiseScale);
    f = f*f*(3.0-2.0*f);
    float n = pf.x + pf.y*57.0 + 113.0*pf.z;

    float h0 = hash(n+0.0);
    float h1 = hash(n+1.0);
    float h57 = hash(n+57.0);
    float h58 = hash(n+58.0);
    float h113 = hash(n+113.0);
    float h114 = hash(n+114.0);
    float h170 = hash(n+170.0);
    float h171 = hash(n+171.0);
    float height = lerp(lerp(lerp(h0,h1,f.x),
                   lerp(h57,h58,f.x),f.y),
               lerp(lerp(h113,h114,f.x),
                   lerp(h170,h171,f.x),f.y),f.z);
    return  clamp(height,-1.0,1.0);
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
// - for some position, returns the distance to the nearest DF object
// ------------------------
vec2 map( in vec3 pos )
{
    vec2 res = opU(vec2(sdCylinder(pos-vec3(0.0,-.25,0.0),vec2(1.9,0.25)),1.0),vec2(sdCylinder(pos,vec2(1.8,0.05)),1.0));
    float pit1 = opS(res.x,sdCylinder(pos, vec2(1.7,0.1)));
    float pit2 = opS(pit1,sdCylinder(pos, vec2(1.6,0.2)));
    float pit3 = opS(pit2,sdCylinder(pos, vec2(1.5,0.3)));
    float pit4 = opS(pit3,sdTorus(pos,vec2(1.0,.5)));
    res = opU( vec2( pit4,1.0),
	                vec2( sdSphere(pos-vec3(-1.0,0.15, 0.0), 0.25 ), 2.0 ) );
    res = opU( res,
               vec2( udRoundBox(pos-vec3(1.0,0.5,0.0),vec3(0.25,0.15,0.25), .05), 3.0) );
    res = opU( res,
               vec2(sdSphere(pos-vec3(0.0,-0.05, 0.0), 0.25 ), 4.0) );
    return res;
}

// ------------------------
// Rendery Bits
// ------------------------
const float PI = 3.14159265;
const vec3 dirLight = normalize(vec3(-5.0,15.0,-5.0)); //this is TO light, kind of
const vec3 ptLight = vec3(2.5,5.0,2.5);
const float tmin = .001;
const float tmax = 10.0;
const float staticT = 0.1;
const int iterations = 100;
const float precis = 0.002;
const float eps = 0.05;
const float k = 16.0; //ao, soft shadows
const vec3 ambient = vec3(0.5,0.5,0.85);
const float ambK = 0.011;
const bool sphereJump = true;

vec3 getMaterialColor(in float id)
{
    if(id == 1.0) {return vec3(.85);}
    if(id == 2.0) {return vec3(.85,0.0,0.0);}
    if(id == 3.0) {return vec3(0.0,.85,0.0);}
    if(id == 4.0) {return vec3(0.0,0.0,.85);}
    else return vec3(0.25);
}

vec3 getBump(in float id, in vec2 uv)
{
    //if(id == 1.0) {return texture2D(iChannel0,uv).xyz;}
    //else if(id == 3.0) {return texture2D(iChannel1,uv).xyz;}
    return vec3(0.0);
}

vec3 castRay( in vec3 ro, in vec3 rd )
{
    float t = tmin;
    float m = -1.0;
    int blah = 0;

    //"Geometry" loop
    for(int i=0; i<iterations; i++)
    {
	    vec2 res = map( ro+rd*t );
        if( res.x<precis || t>tmax ) break;
        if(!sphereJump) {t += staticT;} //constant t increments
        else {t += res.x;} //jump by nearest distance
	    m = res.y;
        blah++;
    }

    //"Terrain" loop - has to be naive, no sphere jumping!
    float t2 = tmin;
    for(int i=0; i<iterations; i++)
    {
        vec3 pt = ro+rd*t2;
        if(pt.y < noise(pt)) {break;}
        else t2 += staticT;
    }

    if(t2 < t) { t = t2; m = 999.0;}
    if( t>tmax ) m=-1.0;
    return vec3( t, m, float(blah));
}

vec3 getNormal(in vec3 p, in float m )
{
    vec2 hxe = map(p + vec3(eps,0,0));
    vec2 hye = map(p + vec3(0,eps,0));
    vec2 hze = map(p + vec3(0,0,eps));

    vec3 n = (1.0/eps) * (map(p).x - vec3(hxe.x,hye.x,hze.x));

    vec2 uv = vec2(0.5 + atan(n.z,n.x)/(2.0 * PI),0.5 - asin(n.y)/PI);
    vec3 alt = getBump(m,uv);

    //return normalize(-n);
    return -normalize(n + alt);
}


float getShadow(in vec3 ro, in vec3 rd) { //ro - point, rd - light direction
    vec3 ray = castRay(ro,rd);
    float t = ray.x;
    float m = ray.y;
    if(m>-0.5) {
        return 0.0;
    }
    return 1.0;
}

float getSoft(in vec3 ro, in vec3 rd) {
    float shade = 1.0;
    float t = tmin;
    for(int i = 0; i < 50; i++) {
        float dist = map(ro + t*rd).x;
        if(dist < eps) {return 0.0;}
        shade = min(shade, k*dist/t);
        if(!sphereJump) {t += staticT;} //constant t increments
        else {t += dist;} //jump by nearest distance
    }
    return clamp(shade, 0.0, 1.0);
}

float getLighting(in vec3 n, in vec3 p) {
    vec3 ptDir = normalize(ptLight - p);
    return dot(n,ptDir) * getSoft(p,ptDir);
}

float getAO(in vec3 p, in float m) {
    //separately calculate the 5 'feelers'?
    //totally should have done this in a loop :(
    vec3 n = getNormal(p,m);
    float f = 0.5 * (eps - map(p + n*eps).x);
    f = f + 0.25 * (2.0*eps - map(p + n*2.0*eps).x);
    f = f + 0.125 * (3.0*eps - map(p + n*3.0*eps).x);
    f = f + 0.0625 * (4.0*eps - map(p + n*4.0*eps).x);
    f = f + 0.03125 * (5.0*eps - map(p + n*5.0*eps).x);
    return clamp(1.0 - k*f,0.0,1.0);
}

vec3 render(in vec3 ro, in vec3 rd, in vec2 coord) {
    //return rd;  // camera ray direction debug view
    //vec3 col = texture2D( iChannel0, coord ).xyz;//vec3(0.1, 0.1, 0.1);
    vec3 col = ambient;
    vec3 res = castRay(ro,rd);
    float t = res.x;
	float m = res.y;
    vec3 p = ro + t*rd;
    p = p - eps*rd;

//toggle for distance-from-camera view
#if 0
    return clamp(vec3(1.0-t/tmax,1.0-t/tmax,0.0),0.0,1.0);
#endif

//toggle for iterations
#if 0
    return clamp(vec3(1.0 - float(res.z)/float(iterations)),0.0,1.0);
#endif

    if( m>-0.5) {
        col = clamp(getMaterialColor(m)
                    * getLighting(getNormal(p,m),p)
                    * getAO(p,m)
                    ,0.0,1.0)
                    + ambK * ambient;
        //col = clamp((ambK*ambient + getMaterialColor(m)*(1.0 - spint) + specular*spint)
        //            * getLighting(getNormal(p),p)
        //            * getAO(p)
        //            ,0.0,1.0);
    }

//toggle for normals render
#if 0
    return vec3( abs(getNormal(p,m)));
#endif

//toggle for AO render
#if 0
    return clamp(vec3(0.85)*getAO(p,m),0.0,1.0);
#endif

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
    vec3 col = render(ro, rd, p);

    col = pow(col, vec3(0.4545));

    fragColor = vec4(col, 1.0);
}
