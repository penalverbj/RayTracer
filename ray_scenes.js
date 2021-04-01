// routines for creating a ray tracing scene
//Jose Penalver Bartolome

let backgroundColor = [1,1,1]; //global array of rgb values for background color
let fov; //global field of view
let eye = {
    x: 0,
    y: 0,
    z: 0,
  }; //global array of xyz position for eye
let u;
let v;
let w;
let ambientLight = [1,1,1]; //global array of rgb values for ambient light
let pointLights = []; //array of point light objects
let spheres = []; //array of sphere objects
let disks = []; //array of disk objects
let areaLights = []; //array of areaLights
let sampleLevel;
let jitter;

// NEW COMMANDS FOR PART B

// create a new disk
function new_disk (x2, y2, z2, radius2, nx2, ny2, nz2, dr2, dg2, db2, k_ambient2, k_specular2, specular_pow2) {
  let d = {
    x: x2,
    y: y2,
    z: z2,
    radius: radius2,
    nx: nx2,
    ny: ny2,
    nz: nz2,
    dr: dr2,
    dg: dg2,
    db: db2,
    k_ambient: k_ambient2,
    k_specular: k_specular2,
    specular_pow: specular_pow2
  };
  
  disks.push(d);
}

// create a new area light source
function area_light (r2, g2, b2, x2, y2, z2, ux2, uy2, uz2, vx2, vy2, vz2) {
  let al = {
    r: r2,
    g: g2,
    b: b2,
    x: x2,
    y: y2,
    z: z2,
    ux: ux2,
    uy: uy2,
    uz: uz2,
    vx: vx2,
    vy: vy2,
    vz: vz2
  };
  
  areaLights.push(al);
}

function set_sample_level (num) {
  sampleLevel = num;
}

function jitter_on() {
  jitter = true;
}

function jitter_off() {
  jitter = false;
}


// OLD COMMANDS FROM PART A (some of which you will still need to modify)


// clear out all scene contents
function reset_scene() {
  backgroundColor = [1,1,1];
  fov = 0;
  eye = {
    x: 0,
    y: 0,
    z: 0,
  };
  let u;
  let v;
  let w;
  ambientLight = [1,1,1];
  pointLights = [];
  spheres = [];
  areaLights = [];
  disks = [];
}

// create a new point light source
function new_light (rp, gp, bp, xp, yp, zp) {
  //creates a pointLight object with necessary variables
  let pl = {
    x: xp,
    y: yp,
    z: zp,
    r: rp,
    g: gp,
    b: bp
  };
  
  pointLights.push(pl);
}

// set value of ambient light source
function ambient_light (r, g, b) {
  ambientLight = [r,g,b];
}

// set the background color for the scene
function set_background (r, g, b) {
  backgroundColor = [r,g,b];
}

// set the field of view
function set_fov (theta) {
  fov = radians(theta);
}

// set the position of the virtual camera/eye
function set_eye_position (x, y, z) {
  eye.x = x;
  eye.y = y;
  eye.z = z;
}

// set the virtual camera's viewing direction
function set_uvw(x1,y1, z1, x2, y2, z2, x3, y3, z3) {
  u = createVector(x1,y1,z1);
  u.normalize();
  v = createVector(x2,y2,z2);
  v.normalize();
  w = createVector(x3,y3,z3);
  w.normalize();
}

// create a new sphere
function new_sphere (x2, y2, z2, radius2, dr2, dg2, db2, k_ambient2, k_specular2, specular_pow2) {
  let s = {
    x: x2,
    y: y2,
    z: z2,
    radius: radius2,
    dr: dr2,
    dg: dg2,
    db: db2,
    k_ambient: k_ambient2,
    k_specular: k_specular2,
    specular_pow: specular_pow2
  };
  
  spheres.push(s);
}

// create an eye ray based on the current pixel's position
function eye_ray_uvw (i, j) {
  //calculates necesary scalars
  let d = 1 / tan(fov/2);
  let u2 = -1 + (2*i)/width;
  let v2 = -1 + (2*(height - j))/height;
  //direction = -dw + uu + vv
  let direction = p5.Vector.add(p5.Vector.add(p5.Vector.mult(w, -d), p5.Vector.mult(u, u2)), p5.Vector.mult(v, v2)); //3d vector
  direction.set(direction.x, direction.y, direction.z); //inverts y so images show up right
  
  let ray = {
   origin: [eye.x, eye.y, eye.z], 
   d: direction
  };
  return ray;
}

//finding the "t" intersections between a ray and a sphere using the quadratic formula
//return -1 if there is no intersection
//return minimum t value intersection with color data and normal vector otherwise
function sphereInter(ray, sphere, subx, suby) { 
  //a = dx^2 + dy^2 + dz^2
  let a = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
  
  //b = 2(x0dx + y0dy + z0dz -cxdx - cydy - czdz)
  let b = ray.origin[0] * ray.d.x + ray.origin[1] * ray.d.y + ray.origin[2] * ray.d.z -
          sphere.x * ray.d.x - sphere.y * ray.d.y - sphere.z * ray.d.z;
  b *= 2;
  
  //c = x0^2 + y0^2 + z0^2 + cx^2 + cy^2 + cz^2 - 2x0cx - 2y0cy - 2z0cz - r^2
  let c = ray.origin[0] * ray.origin[0] + ray.origin[1] * ray.origin[1] + ray.origin[2] * ray.origin[2] +
          sphere.x * sphere.x + sphere.y * sphere.y + sphere.z * sphere.z - 
          2*sphere.x*ray.origin[0] - 2*sphere.y*ray.origin[1] - 2*sphere.z*ray.origin[2] - sphere.radius * sphere.radius;
            
  if ((b*b - 4*a*c) < 0) {
   return -1; 
  } 
  else {
    //calculates both t values
    let t1 = (-b + Math.sqrt(b*b - 4*a*c)) / (2*a);
    let t2 = (-b - Math.sqrt(b*b - 4*a*c)) / (2*a);
    let tin;
    if (t1 < 0 && t2 < 0) {
      return -1;
    }
    //gets smallest t
    else {
      if (t1 < 0) {
        tin = t2;
      }
      if (t2 < 0) {
        tin = t1;
      }
      if (t1 >= 0 && t2 >= 0) {
        tin = min(t1, t2); 
      }
      //calculate normal vector
      let x = (ray.origin[0] + tin*ray.d.x) - sphere.x;
      let y = (ray.origin[1] + tin*ray.d.y) - sphere.y;
      let z = (ray.origin[2] + tin*ray.d.z) - sphere.z;
      let n = createVector(x,y,z);
      n.normalize();
      
      let ls = []; //vectors to point lights from point on the surface
      for (let i = 0; i < pointLights.length; i++) {
       let f = pointLights[i].x - (ray.origin[0] + tin*ray.d.x);
       let s = pointLights[i].y - (ray.origin[1] + tin*ray.d.y);
       let t = pointLights[i].z - (ray.origin[2] + tin*ray.d.z);
       let l = createVector(f,s,t);
       let o = {
         origin: [ray.origin[0] + tin*ray.d.x + 0.0001*n.x, 
                   ray.origin[1] + tin*ray.d.y + 0.0001*n.y, 
                   ray.origin[2] + tin*ray.d.z + 0.0001*n.z],
         d: l,
       };
       ls.push(o); //vector at ls[i] corresponds to pointLights[i]
      }
      //for jitter
      let r = 0;
      if (jitter == true) {
        r = Math.random() - 0.5;
      }
      let al = []; //vectors to point lights from point on the surface
      for (let i = 0; i < areaLights.length; i++) {
        //soft shadow calculations
           U = createVector(areaLights[i].ux, areaLights[i].uy, areaLights[i].uz);
           V = createVector(areaLights[i].vx, areaLights[i].vy, areaLights[i].vz);
           subx2 = ((subx+1+r) / (sampleLevel+1)) * 2 - 1;
           suby2 = ((suby+1+r) / (sampleLevel+1)) * 2 - 1;
           U = p5.Vector.mult(U, subx2);
           V = p5.Vector.mult(V, suby2);
           let f = areaLights[i].x - (ray.origin[0] + tin*ray.d.x) + U.x + V.x; //+ U.x + V.x
           let s = areaLights[i].y - (ray.origin[1] + tin*ray.d.y) + U.y + V.y; //U.y...
           let t = areaLights[i].z - (ray.origin[2] + tin*ray.d.z) + U.z + V.z; 
           let l = createVector(f,s,t);
         let o = {
           origin: [ray.origin[0] + tin*ray.d.x + 0.0001*n.x, 
                   ray.origin[1] + tin*ray.d.y + 0.0001*n.y, 
                   ray.origin[2] + tin*ray.d.z + 0.0001*n.z],
           d: l,
         };
         al.push(o); //vector at al[i] corresponds to areaLights[i]
      }
      
      //return t value, color data, normal vector, light vectors
      //type 0 = sphere intersection
      let intersection = {
        t: tin,
        type: 0,
        normal: n,
        lights: ls,
        Alights: al,
        dr: sphere.dr,
        dg: sphere.dg,
        db: sphere.db,
        k_ambient: sphere.k_ambient,
        k_specular: sphere.k_specular,
        specular_pow: sphere.specular_pow
      };
      
      return intersection;
    }
  }
}

function diskInter(ray, disk, subx, suby) {
 //first check if it hits the plane
 //n = disks normal
 let tin;
 //might have to negate d
 let d = -(disk.nx*disk.x + disk.ny*disk.y + disk.nz*disk.z); 
 let num = -d - (disk.nx*ray.origin[0] + disk.ny*ray.origin[1] + disk.nz*ray.origin[2]);
 let denom = ray.d.x*disk.nx + ray.d.y*disk.ny + ray.d.z*disk.nz;
 
 //intersection found
 if (denom != 0) {
   tin = num / denom;
   if (tin > 0) { //intersection is infront of eye
     //intersection point
     let p = [ray.origin[0] + tin*ray.d.x, ray.origin[1] + tin*ray.d.y, ray.origin[2] + tin*ray.d.z]; 
     let dist = sqrt(sq(disk.x - p[0]) + sq(disk.y - p[1]) + sq(disk.z - p[2]));
     //point is inside the disk
     if (dist <= disk.radius) {
       //POINT light vectors
       let ls = []; //vectors to point lights from point on the surface
         for (let i = 0; i < pointLights.length; i++) {
           let f = pointLights[i].x - (p[0]);
           let s = pointLights[i].y - (p[1]);
           let t = pointLights[i].z - (p[2]);
           let l = createVector(f,s,t);
           let o = {
             origin: [p[0] + 0.0001 * disk.nx, p[1] + 0.0001 * disk.ny, p[2] + 0.0001 * disk.nz],
             d: l
           };
           ls.push(o); //vector at ls[i] corresponds to pointLights[i]
          }
          
          //for jitter
          let r = 0;
          if (jitter == true) {
            r = Math.random() - 0.5;
          }
         //AREA light vectors
         let al = []; //vectors to point lights from point on the surface
         for (let i = 0; i < areaLights.length; i++) {
           //soft shadow calculations
           U = createVector(areaLights[i].ux, areaLights[i].uy, areaLights[i].uz);
           V = createVector(areaLights[i].vx, areaLights[i].vy, areaLights[i].vz);
           subx2 = ((subx+1+r) / (sampleLevel+1)) * 2 - 1;
           suby2 = ((suby+1+r) / (sampleLevel+1)) * 2 - 1;
           U = p5.Vector.mult(U, subx2);
           V = p5.Vector.mult(V, suby2);
           let f = areaLights[i].x - (p[0]) + U.x + V.x; //+ U.x + V.x
           let s = areaLights[i].y - (p[1]) + U.y + V.y; //U.y...
           let t = areaLights[i].z - (p[2]) + U.z + V.z; 
           let l = createVector(f,s,t);
           let o = {
             origin: [p[0] + 0.0001 * disk.nx, p[1] + 0.0001 * disk.ny, p[2] + 0.0001 * disk.nz],
             d: l
           };
           al.push(o); //ray at ls[i] corresponds to areaLights[i]
          }
          
         //normal vector
         let n = createVector(disk.nx, disk.ny, disk.nz); 
         //type 1 = disk intersection
         let i = {
           t: tin,
           type: 1,
           normal: n,
           lights: ls,
           Alights: al,
           dr: disk.dr,
           dg: disk.dg,
           db: disk.db,
           k_ambient: disk.k_ambient,
           k_specular: disk.k_specular,
           specular_pow: disk.specular_pow
         };
         return i;
     } else {
      return -1; 
     }
   } else {
    return -1; 
 }
 } else {
  return -1; 
 }
}

//returns the interesction with the minimum t-value
function findIntersections(ray, subx, suby) {
  //finding fist intersection between all disks and all spheres
      let minT = Number.POSITIVE_INFINITY;
      let min = null;
      //spheres
      for (let i = 0; i < spheres.length; i++) {
        let inter = sphereInter(ray, spheres[i], subx, suby);
        if (inter != -1) { //intersection found
          if (inter.t < minT) { //keeps track of what the first intersection is
            minT = inter.t;
            min = inter;
          }
        }
      }
      //disks
      for (let i = 0; i < disks.length; i++) {
        let inter = diskInter(ray, disks[i], subx, suby);
        if (inter != -1) { //intersection found
          if (inter.t < minT) { //keeps track of what the first intersection is
            minT = inter.t;
            min = inter;
          }
        }
      }
    return min;
}


// this is the main routine for drawing your ray traced scene
function draw_scene() {
  console.log(sampleLevel);
  noStroke();  // so we don't get a border when we draw a tiny rectangle

  // go through all the pixels in the image
  
  for (let y = 0; y < height; y++) {
    for (let x = 0; x < width; x++) {
      let r = 0;  // placeholders to store the pixel's color
      let g = 0;
      let b = 0;
      //subsampling of pixel (x,y)
      for (let sy = 0; sy < sampleLevel; sy++) {
        for (let sx = 0; sx < sampleLevel; sx++) {
          let subx = x + ((sx + 1) / (sampleLevel + 1)) - 0.5;
          let suby = y + ((sy + 1) / (sampleLevel + 1)) - 0.5;
          // create eye ray
          let ray = eye_ray_uvw (subx, suby);
          // maybe print debug information
          debug_flag = 0;
          //if (x == width / 2 && y == height / 2) { debug_flag = 1;  }  // un-comment to debug center pixel
          
          //if (debug_flag) {
          //  console.log ("debug at: " + x + " " + y);
          //}
         
          //Minimum Intersection from eye
          minInter = findIntersections(ray, sx, sy);
          
          if(minInter != null) { //there is an intersection for pixel (i,j)
            r += minInter.k_ambient * minInter.dr * ambientLight[0];
            g += minInter.k_ambient * minInter.dg * ambientLight[1];
            b += minInter.k_ambient * minInter.db * ambientLight[2];  
            
            //loops through all POINT lights updating rgb
            for (let i = 0; i < minInter.lights.length; i++) { 
              shadowT = findIntersections(minInter.lights[i], sx, sy);
              if(shadowT != null) { //shadow needs to be cast
                if (shadowT.t > 0 && shadowT.t < 1) {
                  r += 0;
                  g += 0;
                  b += 0;
                }
              } else { //shadow not casted
                 r += minInter.dr * pointLights[i].r * max(0, p5.Vector.dot(minInter.normal, minInter.lights[i].d.normalize()));
                 g += minInter.dg * pointLights[i].g * max(0, p5.Vector.dot(minInter.normal, minInter.lights[i].d.normalize()));
                 b += minInter.db * pointLights[i].b * max(0, p5.Vector.dot(minInter.normal, minInter.lights[i].d.normalize()));
              }
            }
            //loops through all AREA lights updating rgb
            for (let i = 0; i < minInter.Alights.length; i++) { 
              shadowT = findIntersections(minInter.Alights[i], sx, sy);
              if(shadowT != null) { //shadow needs to be cast
                if (shadowT.t > 0 && shadowT.t < 1) {
                  r += 0;
                  g += 0;
                  b += 0;
                }
              } else { //shadow not casted
                 r += minInter.dr * areaLights[i].r * max(0, p5.Vector.dot(minInter.normal, minInter.Alights[i].d.normalize()));
                 g += minInter.dg * areaLights[i].g * max(0, p5.Vector.dot(minInter.normal, minInter.Alights[i].d.normalize()));
                 b += minInter.db * areaLights[i].b * max(0, p5.Vector.dot(minInter.normal, minInter.Alights[i].d.normalize()));
              }
            }
          } else { //no intersection so just background color
            r += backgroundColor[0];
            g += backgroundColor[1];
            b += backgroundColor[2]; 
          }
        }
      }
      //averages rgb values for anti-aliasing  
      r = r / sq(sampleLevel);
      g = g / sq(sampleLevel);
      b = b / sq(sampleLevel);
      //set the pixel color, converting values from [0,1] into [0,255]
      fill (255 * r, 255 * g, 255 * b);
      
      rect (x, y, 1, 1);   // make a little rectangle to fill in the pixel
    }
  }
  
}
