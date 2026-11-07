import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title='Plane Stress Explorer', layout='wide')
st.title('Plane Stress Transformation Explorer')
st.caption('Rotate a material element, inspect its stress components, and compare with Mohr’s circle. Units: MPa and degrees.')

with st.sidebar:
    st.header('Original stress state')
    sx = st.number_input('σx (MPa)', value=50.0, step=5.0)
    sy = st.number_input('σy (MPa)', value=-10.0, step=5.0)
    txy = st.number_input('τxy (MPa)', value=40.0, step=5.0)

avg = (sx + sy)/2
half = (sx - sy)/2
radius = np.hypot(half, txy)
principal = np.degrees(0.5*np.arctan2(2*txy, sx-sy)) if radius > 1e-12 else 0.0
shear_angle = principal-45

def transform(theta):
    a=np.radians(2*theta)
    d=half*np.cos(a)+txy*np.sin(a)
    return avg+d, avg-d, -half*np.sin(a)+txy*np.cos(a)

st.subheader('Rotation')
mode=st.radio('Jump to an orientation', ['Manual', 'Maximum principal stress', 'Minimum principal stress', 'Positive maximum shear', 'Negative maximum shear'], horizontal=True)
if mode=='Manual':
    theta=st.slider('θ (degrees, counterclockwise)', -180.0,180.0,0.0,0.5)
else:
    theta={'Maximum principal stress':principal,'Minimum principal stress':principal+90,'Positive maximum shear':shear_angle,'Negative maximum shear':shear_angle+90}[mode]
    st.info(f'θ = {theta:.2f}°')
sxp,syp,tp=transform(theta)
c1,c2,c3=st.columns(3)
c1.metric('σx′', f'{sxp:.2f} MPa')
c2.metric('σy′', f'{syp:.2f} MPa')
c3.metric('τx′y′', f'{tp:.2f} MPa')

fig, (ax,mc)=plt.subplots(1,2,figsize=(12,5))
ang=np.radians(theta)
ex=np.array([np.cos(ang),np.sin(ang)])
ey=np.array([-np.sin(ang),np.cos(ang)])
scale=0.9
corners=np.array([-ex-ey,ex-ey,ex+ey,-ex+ey,-ex-ey])*scale
ax.plot(corners[:,0],corners[:,1],color='steelblue',lw=2)
ax.fill(corners[:,0],corners[:,1],color='steelblue',alpha=.12)
maxstress=max(abs(sxp),abs(syp),abs(tp),1)

def draw_face(normal,tangent,normalstress,shearstress):
    center=normal*scale
    # Positive normal stress pulls outward; negative pushes inward.
    for val, direction, color, offset in [(normalstress,normal,'darkorange',0),(shearstress,tangent,'teal',.25)]:
        if abs(val)<1e-9: continue
        start=center + tangent*offset
        length=0.72*abs(val)/maxstress
        direction=direction*np.sign(val)
        ax.arrow(start[0],start[1],direction[0]*length,direction[1]*length,
                 head_width=.075,head_length=.11,length_includes_head=True,
                 color=color,lw=2)
for n,t,ns,ts in [(ex,ey,sxp,tp),(-ex,-ey,sxp,tp),(ey,-ex,syp,-tp),(-ey,ex,syp,-tp)]:
    draw_face(n,t,ns,ts)
ax.arrow(0,0,1.55*ex[0],1.55*ex[1],head_width=.06,color='gray',alpha=.65)
ax.arrow(0,0,1.55*ey[0],1.55*ey[1],head_width=.06,color='gray',alpha=.65)
ax.text(*(1.7*ex),'x′',fontsize=12)
ax.text(*(1.7*ey),'y′',fontsize=12)
ax.axhline(0,color='gray',alpha=.2);ax.axvline(0,color='gray',alpha=.2)
ax.set_aspect('equal');ax.set_xlim(-2.5,2.5);ax.set_ylim(-2.5,2.5)
ax.set_title('Rotated element (orange: normal; teal: shear)')
ax.set_xlabel('Global x');ax.set_ylabel('Global y')

circle_t=np.linspace(0,2*np.pi,500)
mc.plot(avg+radius*np.cos(circle_t),radius*np.sin(circle_t),label='Mohr circle')
mc.scatter([sxp,syp],[tp,-tp],c=['darkorange','teal'],zorder=5)
mc.annotate("(σx′, τx′y′)",(sxp,tp),xytext=(8,8),textcoords='offset points')
mc.annotate("(σy′, −τx′y′)",(syp,-tp),xytext=(8,-15),textcoords='offset points')
mc.scatter([avg+radius,avg-radius],[0,0],color='crimson',zorder=6)
mc.axhline(0,color='gray',lw=.7);mc.axvline(0,color='gray',lw=.7)
mc.set_aspect('equal',adjustable='datalim');mc.grid(alpha=.2)
mc.set_xlabel('Normal stress σ (MPa)');mc.set_ylabel('Shear stress τ (MPa)')
mc.set_title('Mohr circle (τ positive upward)')
st.pyplot(fig)
plt.close(fig)
st.caption('Arrow lengths are proportional to stress magnitude. Opposite faces have opposite traction vectors. Mohr-circle angle changes by twice the physical rotation angle, with opposite sense in this convention.')

st.subheader('Analytical results')
a,b,c,d=st.columns(4)
a.metric('σ₁',f'{avg+radius:.2f} MPa')
b.metric('σ₂',f'{avg-radius:.2f} MPa')
c.metric('|τ|max, in-plane',f'{radius:.2f} MPa')
d.metric('Normal stress on max-shear planes',f'{avg:.2f} MPa')
st.write(f'Principal planes: θp = {principal:.2f}° and {principal+90:.2f}° (modulo 180°).')
st.write(f'Max-shear planes: θs = {shear_angle:.2f}° and {shear_angle+90:.2f}° (modulo 180°).')
st.info('In plane stress, the absolute 3D maximum shear may differ: τmax,3D = [max(σ₁, σ₂, 0) − min(σ₁, σ₂, 0)] / 2.')
with st.expander('Equations and derivation',expanded=True):
    st.latex(r'\sigma_{x\prime}=\frac{\sigma_x+\sigma_y}{2}+\frac{\sigma_x-\sigma_y}{2}\cos 2\theta+\tau_{xy}\sin 2\theta')
    st.latex(r'\sigma_{y\prime}=\frac{\sigma_x+\sigma_y}{2}-\frac{\sigma_x-\sigma_y}{2}\cos 2\theta-\tau_{xy}\sin 2\theta')
    st.latex(r'\tau_{x\prime y\prime}=-\frac{\sigma_x-\sigma_y}{2}\sin 2\theta+\tau_{xy}\cos 2\theta')
    st.latex(r'\sigma_{1,2}=\frac{\sigma_x+\sigma_y}{2}\pm\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}')
    st.latex(r'\tau_{\max,\mathrm{in\ plane}}=\sqrt{\left(\frac{\sigma_x-\sigma_y}{2}\right)^2+\tau_{xy}^2}')
    st.latex(r'2\theta_p=\operatorname{atan2}(2\tau_{xy},\sigma_x-\sigma_y)')
    st.latex(r'\theta_s=\theta_p\pm45^\circ')
    st.write('Principal orientation comes from setting transformed shear to zero. Maximum shear orientation comes from setting the derivative of transformed shear with respect to θ to zero. atan2 resolves the correct quadrant, unlike tan⁻¹ of a simple ratio.')
    if radius<1e-12:
        st.warning('Equal biaxial stress: every in-plane orientation is principal and shear is zero; principal angles are not unique.')
