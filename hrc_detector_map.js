//JavaScript code for simulation of X-ray Laue backscattering

var version = "0.1";

// dimensions of the canvas object
var scaleX=1200;
var scaleY=400;

var X0=scaleX/2;
var Y0=scaleY/2;
var X0_ofst=0;
var Y0_ofst=0;


//parameters for the appearance of the simulation
const decimal_digit=1000;     // decimal digit for UBmatrix
const radius=5;       // radius of circles showing refletions in the simulation.
const radius_tgt=8;     //// radius of a circle showing a target refletions in the simulation.
const txt_ofst1=radius+10;   //offset along Y direction for indices shown near each reflection.
const txt_ofst2=3;   //offset along X direction for detector number shown bottom.
const fundamental_color="rgb(0, 0, 250)";
const DetMapBGColor="rgb(220, 220, 220)";
const gridcolor="rgb(250, 100, 0)"
const ref_linewidth=1;

//variables for calculating Laue diffraction patterns.-------------------
var u = new Array(3); // indices, pallarel to the incident beam
var v = new Array(3); // indices, another direction in the horizontal plane including the incidnet beam

var ux = new Array(3);
var vx = new Array(3);

var Rot0 = new Array(3);
var Rot1 = new Array(3);
var Rot2 = new Array(3);
var Rot =[Rot0, Rot1, Rot2];    // 3x3 rotation matrix

var a_unit = new Array(3);  // unit vector of primitive translation vectors
var b_unit = new Array(3);
var c_unit = new Array(3);

var a_star = new Array(3);  // reciprocal lattice vectors
var b_star = new Array(3);
var c_star = new Array(3);

var as_len;
var bs_len;
var cs_len;

var RefCon = '';

var Hmax;
var Kmax;
var Lmax;

var lambda_min=0.4;

var phih;
var phiv;
var lambda;             // wavelength 

var Omega=0;

//parameters regarding the detector banks
var Lsd = 40;   // Distance between the sample and detector (mm)
var DetHeight = 80; //height of the detector (mm)
var HD = 20;    // height of center of PSD from incident beam (mm)
var LD = 2800;  // length of PSD (mm)

//variables for 3D orientation viewer
const arrow_scale = 120;        //arrows for a*, b* and c*: convert A-1 to pixel.
const arrow_HeadLen = 20;       //lengths of arrowheads (pixel)
const arrow_HeadWidth = 10;     //widths of arrowheads (pixel)
const scale3D = 5;   // convert mm to pixel.
const DetBankThickness = 50; //pixel

//variable for loading observed Laue image.
var imageLoaded=false;
var imageURL;
var image = new Image();

function init_draw(){
    document.getElementById("verNum").innerHTML=version;
    document.getElementById("verNum2").innerHTML=version;
    draw();
}

function draw() {
    set_Lattice();
    set_ReflectionCondition();
    lambda_adjust_and_draw();
    draw_OriViewer();

}

function set_RefCon_and_draw(){
    set_ReflectionCondition();
    draw_DetMap();
}

function rot_and_draw(rot_ax_dir) {
    rot_Lattice(rot_ax_dir);
    draw_DetMap();
    draw_OriViewer();
}

function lambda_adjust_and_draw(){
    document.getElementById("lambda_min_disp").value = document.getElementById("lambda_min").value;
    lambda_min = Number(document.getElementById("lambda_min").value);
    draw_DetMap();
}

function set_Origin_and_draw(){
    X0_ofst = Number(document.getElementById("X0_ofst").value);
    Y0_ofst = -Number(document.getElementById("Y0_ofst").value);
    Lsd = Number(document.getElementById("Lsd").value);
    DetHeight = Number(document.getElementById("DetHeight").value);
    draw_DetMap();
}

function omegaRot_and_draw(){
    let targetOmega=Number(document.getElementById("targetOmega").value);
    const deltaOmega = (targetOmega-Omega)/180.0*Math.PI;
    xyz_rotation(2,-deltaOmega);    //"2" means rotation about the z axis. a minus sign is necessary because directions of omega- and z-rotations are opposite to each other.
    draw_DetMap();
    draw_OriViewer();
    Omega=targetOmega;
    document.getElementById("currentOmega_disp").innerHTML=String(Omega);
}

function set_Lattice(){

    //input parameters: lattice constants and sample orientation)
    let a = Number(document.getElementById('a').value);
    let b = Number(document.getElementById('b').value);
    let c = Number(document.getElementById('c').value);
    let alpha = Number(document.getElementById('alpha').value)/180.0*Math.PI;   // in radian
    let beta  = Number(document.getElementById('beta').value)/180.0*Math.PI;    // in radian
    let gamma = Number(document.getElementById('gamma').value)/180.0*Math.PI;   // in radian
    u[0] = Number(document.getElementById('u1').value);
    u[1] = Number(document.getElementById('u2').value);
    u[2] = Number(document.getElementById('u3').value);
    v[0] = Number(document.getElementById('v1').value);
    v[1] = Number(document.getElementById('v2').value);
    v[2] = Number(document.getElementById('v3').value);

    // calculation
    let DD = (Math.cos(alpha)-Math.cos(gamma)*Math.cos(beta))/Math.sin(gamma);
    let PP = Math.sqrt(Math.sin(beta)-DD**2.0);

    ux[0] = 2.0*Math.PI*u[0]/a;
    ux[1] = 2.0*Math.PI*(-u[0]/a/Math.tan(gamma)+u[1]/b/Math.sin(gamma));
    ux[2] = 2.0*Math.PI*(u[0]/a*(DD/Math.tan(gamma)-Math.cos(beta))-u[1]/b*DD/Math.sin(gamma)+u[2]/c)/PP;
    vx[0] = 2.0*Math.PI*v[0]/a;
    vx[1] = 2.0*Math.PI*(-v[0]/a/Math.tan(gamma)+v[1]/b/Math.sin(gamma));
    vx[2] = 2.0*Math.PI*(v[0]/a*(DD/Math.tan(gamma)-Math.cos(beta))-v[1]/b*DD/Math.sin(gamma)+v[2]/c)/PP;

    let uy2uz2 = ux[1]**2.0+ux[2]**2.0;    
    let Uabs = Math.sqrt(ux[0]**2.0+uy2uz2);
    let Rvy;
    let Rvz;
    if(uy2uz2==0){
        Rvy=vx[1];
        Rvz=vx[2];
    }
    else{
        Rvy =(-vx[0]*ux[1]+(vx[1]*(ux[0]*ux[1]**2.0+Uabs*ux[2]**2.0)+vx[2]*ux[1]*ux[2]*(ux[0]-Uabs))/uy2uz2)/Uabs;
        Rvz =(-vx[0]*ux[2]+(vx[2]*(ux[0]*ux[2]**2.0+Uabs*ux[1]**2.0)+vx[1]*ux[2]*ux[1]*(ux[0]-Uabs))/uy2uz2)/Uabs;
    }
    
    let cosphi=Rvy/Math.sqrt(Rvy**2.0+Rvz**2.0);
    let sinphi=Rvz/Math.sqrt(Rvy**2.0+Rvz**2.0);

    Rot[0][0]= ux[0]/Uabs;
    Rot[0][1]= ux[1]/Uabs;
    Rot[0][2]= ux[2]/Uabs;
    Rot[1][0]= -(ux[1]*cosphi+ux[2]*sinphi)/Uabs;
    Rot[1][1]=(ux[2]*(ux[2]*cosphi-ux[1]*sinphi)+ux[0]*ux[1]*(ux[1]*cosphi+ux[2]*sinphi)/Uabs)/uy2uz2;
    Rot[1][2]=(ux[1]*(ux[1]*sinphi-ux[2]*cosphi)+ux[0]*ux[2]*(ux[2]*sinphi+ux[1]*cosphi)/Uabs)/uy2uz2;
    Rot[2][0]=(ux[1]*sinphi-ux[2]*cosphi)/Uabs;
    Rot[2][1]=(-ux[2]*(ux[1]*cosphi+ux[2]*sinphi)+ux[0]*ux[1]*(ux[2]*cosphi-ux[1]*sinphi)/Uabs)/uy2uz2;
    Rot[2][2]=(ux[1]*(ux[1]*cosphi+ux[2]*sinphi)+ux[0]*ux[2]*(ux[2]*cosphi-ux[1]*sinphi)/Uabs)/uy2uz2;

    for (let i=0;i<3;i++){
        a_unit[i]= Rot[i][0];
        b_unit[i]= Rot[i][0]*Math.cos(gamma)+Rot[i][1]*Math.sin(gamma);
        c_unit[i]= Rot[i][0]*Math.cos(beta)+Rot[i][1]*DD+Rot[i][2]*PP;
    } 

    // output parameters: 3 reciprocal lattice vectors, a*, b*, and c*
    for (let i=0;i<3;i++){
        a_star[i]= 2.0*Math.PI/a/PP/Math.sin(gamma)*(b_unit[(i+1)%3]*c_unit[(i+2)%3]-b_unit[(i+2)%3]*c_unit[(i+1)%3]);
        b_star[i]= 2.0*Math.PI/b/PP/Math.sin(gamma)*(c_unit[(i+1)%3]*a_unit[(i+2)%3]-c_unit[(i+2)%3]*a_unit[(i+1)%3]);
        c_star[i]= 2.0*Math.PI/c/PP/Math.sin(gamma)*(a_unit[(i+1)%3]*b_unit[(i+2)%3]-a_unit[(i+2)%3]*b_unit[(i+1)%3]);
    }
    
    as_len = Math.sqrt(a_star[0]**2.0+a_star[1]**2.0+a_star[2]**2.0);
    bs_len = Math.sqrt(b_star[0]**2.0+b_star[1]**2.0+b_star[2]**2.0);
    cs_len = Math.sqrt(c_star[0]**2.0+c_star[1]**2.0+c_star[2]**2.0);

}

function set_ReflectionCondition(){
    RefCon = document.getElementById("RefCon").value;
}

function check_ReflectionCondition(RefCon,H,K,L){
    var retstr=false;

    switch(RefCon){
        case 'none':
            retstr=true;
            break;
        case 'H+K=2n':
            if((H+K)%2==0){
                retstr=true;
            }
            break;
        case 'H+L=2n':
            if((H+L)%2==0){
                retstr=true;
            }
            break;
        case 'K+L=2n':
            if((K+L)%2==0){
                retstr=true;
            }
            break;
        case 'H+K+L=2n':
            if((H+K+L)%2==0){
                retstr=true;
            }
            break;
        case 'H,K,L all even or all odd':
            let hklsp = Math.abs(H%2)+Math.abs(K%2)+Math.abs(L%2);
            if(hklsp==0||hklsp==3){
                retstr=true;
            }
            break;
        case '-H+K+L=3n':
            if((-H+K+L)%3==0){
                retstr=true;
            }
            break;
        default:
            retstr=true;
    }
    return retstr;
}


function draw_DetMap(){

    var canvas = document.getElementById('CanvasDetMap');
    canvas.width=scaleX;
    canvas.height=scaleY;

    var context = canvas.getContext('2d');

    //refresh
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.strokeStyle = "rgb(0, 0, 0)";
    context.lineWidth=1;
    

    //set background color
    context.fillStyle = DetMapBGColor;
    context.fillRect(0, 0, canvas.width, canvas.height);

    //show observed Laue pattern
    if(imageLoaded==true){
        context.drawImage(image, 0, 0);
    }

    context.strokeStyle = gridcolor;
    context.beginPath();
    context.moveTo(0,Y0+Y0_ofst);
    context.lineTo(scaleX,Y0+Y0_ofst);
    context.stroke();

    context.strokeStyle = gridcolor;
    context.beginPath();
    context.moveTo(X0+X0_ofst,0);
    context.lineTo(X0+X0_ofst,scaleY);
    context.stroke();

    // color setting for circles indicating reflections
    context.strokeStyle = fundamental_color;
    context.fillStyle = fundamental_color;
    context.lineWidth= ref_linewidth;
    context.font = "10px sans-serif";

    let Qmax = 4.0*Math.PI/lambda_min;
    Hmax = Math.floor(Qmax/as_len);
    Kmax = Math.floor(Qmax/bs_len);
    Lmax = Math.floor(Qmax/cs_len);

    let Ghkl=new Array(3);
    let isTargetHKL=false;
    let showHKL=false;
    for (var H=-Hmax;H<=Hmax;H+=1){
        for (var K=-Kmax;K<=Kmax;K+=1){
            for (var L=-Lmax;L<=Lmax;L+=1){

                if(check_ReflectionCondition(RefCon,H,K,L)==false){
                    // Reflection condition is not satisfied or H=K=L=0.
                }
                else{
                    // funcamental reflections
                    context.strokeStyle = fundamental_color;
                    context.fillStyle = fundamental_color;
                    showHKL=false;
                    drawBraggReflection(context,H,K,L,isTargetHKL,showHKL);

                }
            }
        }
    }

    //draw large circle for the target reflection.
    isTargetHKL=true;
    showHKL=true;
    context.strokeStyle = fundamental_color;
    let Ht=-Number(document.getElementById("Ht").value);
    let Kt=-Number(document.getElementById("Kt").value);
    let Lt=-Number(document.getElementById("Lt").value);
    //minus signs are necessary to convert Q=kf-ki to Q=ki-kf.
    drawBraggReflection(context,Ht,Kt,Lt,isTargetHKL,showHKL);


}

function drawBraggReflection(context1,H1,K1,L1,isTargetHKL1,showHKL1){

    let return_value=false;

    let Ghkl=new Array(3);

    for(let i=0;i<3;i++){
        Ghkl[i]=H1*a_star[i]+K1*b_star[i]+L1*c_star[i];
    }
    if(Ghkl[0]>=0.0){
        // Bragg's law is not satisfied.
    }
    else{  
        let G_sq = Ghkl[0]**2.0+Ghkl[1]**2.0+Ghkl[2]**2.0;
        let ki = -0.5*G_sq/Ghkl[0]; // Ki >0
        lambda = 2.0*Math.PI/ki;    // Angstrome

        if(lambda>lambda_min && G_sq > 0){   // lambda_min=2PI/sqrt(Ei_max/2.072), the case that H=K=L=0 is avoided by the condition of  G_sq > 0.
            let kf=new Array(3);
            kf[0]=Ghkl[0]+ki;
            kf[1]=Ghkl[1];
            kf[2]=Ghkl[2];


            if(kf[0]<0){    //Backscattering condition: the scattered X-ray must go to -x direction, namely opposite to the incident x-ray direction. 
                const mm2pixel=scaleY/DetHeight;
                let PosX=-kf[1]/kf[0]*Lsd*mm2pixel+X0+X0_ofst;
                let PosY=kf[2]/kf[0]*Lsd*mm2pixel+Y0+Y0_ofst;
    
                if(PosX>=0 && PosX<scaleX && PosY >= 0 && PosY <=scaleY){
                    context1.beginPath();
                    if(isTargetHKL1==true){
                        context1.arc(PosX,PosY, radius_tgt, 0, 2 * Math.PI);
                    }
                    else{
                        context1.arc(PosX,PosY, radius, 0, 2 * Math.PI);
                    }
                    context1.stroke();
    
                    if(showHKL1==true){
                        context1.fillText(String(-H1)+String(-K1)+String(-L1), PosX, PosY+txt_ofst1);   
                        //Thus far, the Bragg conditions are calculated assuming that the scattering vector is defined as Q=kf-ki.
                        //However, the correct definition of the scattering vector is Q=ki-kf, which is momentum of the excitation. Thus, "-" signs are necessary to change Q=kf-ki to Q=ki-kf.
                    }
    
                    return_value=true;
                }    
            }
        }
    }
    
    return return_value;

}


function rot_Lattice(rot_ax_dir){
    let angle = 0.0;  // radian unit
    let xyz         // xyz=(0,1,2) for (x, y, z)-axis respectively.
    switch(rot_ax_dir){
        case 'rot_x_plus':
            angle = Number(document.getElementById('rot_x_deg').value)/180.0*Math.PI;
            xyz =0.0;
            break;
        case 'rot_x_minus':
            angle = (-1.0)*Number(document.getElementById('rot_x_deg').value)/180.0*Math.PI;
            xyz =0.0;
            break;
        case 'rot_y_plus':
            angle = Number(document.getElementById('rot_y_deg').value)/180.0*Math.PI;
            xyz =1.0;
            break;
        case 'rot_y_minus':
            angle = (-1.0)*Number(document.getElementById('rot_y_deg').value)/180.0*Math.PI;
            xyz =1.0;
            break;
        case 'rot_z_plus':
            angle = Number(document.getElementById('rot_z_deg').value)/180.0*Math.PI;
            xyz =2.0;
            break;
        case 'rot_z_minus':
            angle = (-1.0)*Number(document.getElementById('rot_z_deg').value)/180.0*Math.PI;
            xyz =2.0;
            break;
        default:
    }
    xyz_rotation(xyz,angle);
}

function xyz_rotation(xyz,angle){
    //xyz : 0=x, 1=y, 2=z
    //angle : rotation angle (radian units)
    let r00;
    let r01;

    r00=a_star[(xyz+1)%3]*Math.cos(angle)-a_star[(xyz+2)%3]*Math.sin(angle);
    r01=a_star[(xyz+1)%3]*Math.sin(angle)+a_star[(xyz+2)%3]*Math.cos(angle);
    a_star[(xyz+1)%3]=r00;
    a_star[(xyz+2)%3]=r01;
    r00=b_star[(xyz+1)%3]*Math.cos(angle)-b_star[(xyz+2)%3]*Math.sin(angle);
    r01=b_star[(xyz+1)%3]*Math.sin(angle)+b_star[(xyz+2)%3]*Math.cos(angle);
    b_star[(xyz+1)%3]=r00;
    b_star[(xyz+2)%3]=r01;
    r00=c_star[(xyz+1)%3]*Math.cos(angle)-c_star[(xyz+2)%3]*Math.sin(angle);
    r01=c_star[(xyz+1)%3]*Math.sin(angle)+c_star[(xyz+2)%3]*Math.cos(angle);
    c_star[(xyz+1)%3]=r00;
    c_star[(xyz+2)%3]=r01;
}

function draw_OriViewer(){
    // サイズを指定
    const width = 800;
    const height = 400;
  
    // レンダラーを作成
    const renderer = new THREE.WebGLRenderer({
      canvas: document.querySelector('#OrientationViewer'),
      antialias: true
    });
    renderer.setPixelRatio(window.devicePixelRatio);
    renderer.setSize(width, height);
    renderer.setClearColor(0xf8f8f8);
  
    // シーンを作成
    const scene = new THREE.Scene();
  
    // カメラを作成
    const camera = new THREE.PerspectiveCamera(30, width / height);
    let cam_theta=Number(document.getElementById("cam_theta").value);
    let cam_phi=Number(document.getElementById("cam_phi").value);
    let cam_len=1200;
    camera.position.set(cam_len*Math.sin(Math.PI/180.0*cam_theta)*Math.sin(Math.PI/180.0*cam_phi), cam_len*Math.cos(Math.PI/180.0*cam_theta), cam_len*Math.sin(Math.PI/180.0*cam_theta)*Math.cos(Math.PI/180.0*cam_phi));
    camera.lookAt(new THREE.Vector3(0, 0, 0));
  
    // detector
    const material1 = new THREE.MeshStandardMaterial({ color: 0xC0C0C0 });  // color of detector bank
    let geometry_det = new THREE.BoxGeometry(DetBankThickness,DetHeight*scale3D,DetHeight/scaleY*scaleX*scale3D);
    let mesh_det = new THREE.Mesh(geometry_det, material1);
    scene.add(mesh_det);
    mesh_det.position.x -= Lsd*scale3D;

    // guide for the incident beam
    const geometry_guide = new THREE.BoxGeometry(1000,50,50);
    const mesh_guide = new THREE.Mesh(geometry_guide, material1);
    scene.add(mesh_guide);
    mesh_guide.position.x -= Lsd*scale3D+500;
    
    //draw a*, b*, c*
    //a*
    var dir = new THREE.Vector3( a_star[0],a_star[2], -a_star[1] );
    var origin = new THREE.Vector3( 0, 0, 0 );
    var arrow_len = dir.length()*arrow_scale;
    var hex = 0xff0000;
    var arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin, arrow_len, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);
  
    //b*
    dir = new THREE.Vector3( b_star[0],b_star[2], -b_star[1] );
    arrow_len = dir.length()*arrow_scale;
    hex = 0x00ff00;
    arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin, arrow_len, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);
 
    //c*
    dir = new THREE.Vector3( c_star[0],c_star[2], -c_star[1] );
    arrow_len = dir.length()*arrow_scale;
    hex = 0x0000ff;
    arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin,arrow_len, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);
  
    //ki*
    dir = new THREE.Vector3( 1,0, 0 );
    arrow_len = dir.length()*arrow_scale;
    hex = 0xff9900;
    var origin2 = new THREE.Vector3( -Lsd*scale3D, 0, 0 );
    arrowHelper = new THREE.ArrowHelper( dir.normalize(), origin2,Lsd*scale3D, hex ,arrow_HeadLen,arrow_HeadWidth);
    scene.add(arrowHelper);

    
    // 平行光源
    const directionalLight = new THREE.DirectionalLight(0xffffff);
    directionalLight.position.set(150, 240, -500);
    scene.add(directionalLight);
  
    const light = new THREE.AmbientLight(0xa0a0a0, 1.0);
    scene.add(light);  
    // ポイント光源
  //  const pointLight = new THREE.PointLight(0xffffff, 2, 1000);
  //  scene.add(pointLight);
  //  const pointLightHelper = new THREE.PointLightHelper(pointLight, 3);
  //  scene.add(pointLightHelper);
  
    renderer.render(scene, camera);
  
  }

function getFile(e){
    let reader = new FileReader();
    reader.readAsDataURL(e[0]);
    reader.onload = function() {
        imageLoaded=true;
        image.src = reader.result;
        image.onload=function() {
            scaleX=image.width;
            scaleY=image.height;
            X0=scaleX/2;
            Y0=scaleY/2;
            draw_DetMap();
            draw_OriViewer();         
        };
    };
}


function removeFile(){
    imageLoaded=false;
    draw_DetMap();         
}



//--------------------------------------

