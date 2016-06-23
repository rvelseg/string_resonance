from pylab import plot, show, subplots
from numpy import arange, pi, ma, tanh, sin, fft, zeros
from matplotlib import animation
from detect_peaks import detect_peaks

__author__ = "R. Velasco-Segura and Pablo L. Rendon"
__version__ = "0.2"
__license__ = "BSD"

# ------------------------------------------------------
# ------- configurable parameters

ampl = 1      # amplitude of the emitted wave

L = 3.0       # resonator length in meters
x_max = 4.0   # domain length in meters

Nx = 50                          # number of grid points
dx = x_max/(Nx - 1.0)            # spatial step in meters
x = arange(0,x_max+dx/2.0,dx)    # spatial domain

t_max = 1.0          # final time in seconds

f_ini = 100.0        # initial frequency in Hz
f_fin = 300.0        # final frequency in Hz

c = 343.0             # speed of sound in m/s

# 343 m/s is typically the speed of sound in air, and not necessarily
# on a string. However, it is not that different from typicall values
# in guitar strings.

mic_pos = x_max/8.0     # microphone position in meters

# Chirp type. The possible values this variable can take are:
#
# "exp" : stands for exponential, which sometimes is also called
# logarithmic chirp.
#
# "lin" : linear chirp.
#
chirp_type = "exp"

# Temporal step. The choice dt=dx/c has some advantages. However, be
# aware that in this case the obstacle as harmonic oscillator only is
# not usable: small values of KoM have no effect, and large values of
# KoM make the numeric system unstable, without an intermediate
# region. When the obstacle is only an harmonic oscillator, with no
# estra mass, dt=0.9*dx/c works fine.
dt = 1*dx/c

# Output type. The possible values this variable can take are:
#
# "file" : The script creates a video file. You probably will have to
# adjust the writer to something that works on your system.
#
# "fly" : A visualization is displayed on the fly, as data is been
# generated.
#
output = "fly"

# ---------------------------------------------------
# ------- calculated and fixed paramenters, ---------
# ------- don't change these parameters.

t_min = 0.0            # initial time in seconds
w = 2.0*pi*f_ini       # current value of angular frequency
phi = 0.0              # phase correction
s_fc = 1.0
n = 0                  # current step
init_sim = 0

cdtdx2 = (c*dt/dx)**2
dt2 = dt*dt

t_axis = arange(t_min,t_max+dt/2.0,dt)    # temporal domain

mic_i = int(round(mic_pos/dx))     # microphone position index
obs_i = int(round(L/dx))           # obstacle position index

# ----------------------------------------------------

y_next = zeros(Nx)           # solution at t
y_now = y_next.copy()        # solution at t-dt
y_prev = y_now.copy()        # solution at t-2*dt

mic = []
mic_t = []

am_min = []
am_max = []
am = []
am_s = []

f_ax = []

fft_spec = []
fft_axis = []
fft_axis_t = []
fft_peaks = []

fig, ax = subplots(3)

def init_y() :

    global y_now, y_prev, x

    for i in range(0,len(x)) :
        y_now[i] = 0
        y_prev[i] = 0

def step_y() :

    global y_next, y_now, y_prev
    global dt2, cdtdx2, dx
    global obs_i

    for i in range(1, len(y_next)-1) :
        y_next[i] = (- y_prev[i] +
                 2 * y_now[i] +
                 cdtdx2 * (y_now[i-1] - 2*y_now[i] + y_now[i+1]) )

    # Obstacle :
    i = obs_i

    # # harmonic oscillator
    # KoM = 1e5  # this is k/m
    # y_next[i] -= (dt2/dx)*KoM*y_now[i]

    # extra mass
    MM = 5.0
    kappa = 4e5
    k = kappa/dx
    y_next[i] += (dt2*(k/MM) - cdtdx2) * (y_now[i-1] - 2*y_now[i] + y_now[i+1])

def boundary_y(t) :

    global y_now, y_prev
    global x, w, phi, ampl

    # absorbing boundary
    y_now[-1] = y_prev[-2]

    # forcing boundary
    y_now[0] = ( ampl *
             sin( w * (t - x[0]) + phi) *
             ( 0.5 + 0.5 * tanh( t * 1715 - 4 ) ) )

def pp(x,t) :
    # This function is used to draw a sine wave of the current emitted
    # frequency, to be compared with y_now.
    global w, phi
    ret = ( ampl *
            sin( w * (t - x) + phi) *
            ( 0.5 + 0.5 * tanh( t * 1715 - 4 ) ) )
    return ret

def get_mic(pos, t) :

    global y_now
    global mic, mic_t

    mic.append(y_now[pos])
    mic_t.append(t)

# static like variables for this function
# norm = 0
# prev_s = 1
# prev_dd = 0
def get_am(pos, t) :

    global am_max, am_min
    global w, dt, mic
    # global s_fc
    # global norm, prev_s, prev_dd
    global f_fin

    ii = len(mic)
    period = 2*pi/w
    i_p = int(period/dt)
    i_p_m = max(0, int(ii - 5*i_p))
    # i_p_m2 = max(0, int(ii - 10*i_p))

    am_min.append(min(mic[i_p_m:]))
    am_max.append(max(mic[i_p_m:]))
    am.append(am_max[-1] - am_min[-1])

    f_ax.append(w/(2*pi))

 #    dd = am[-1] - 2*am[i_p_m] + am[i_p_m2]
 #    print dd

 #    print am[-1], am[i_p_m]
 #    t_adp_start = t_min + (t_max-t_min)/10

 #    if (t > t_adp_start ) :

 #        if norm == 0 :
 #            norm = abs(am[-1] - am[i_p_m])/am[-1]

 #        m = abs(am[-1] - am[i_p_m])/(am[-1]*0.5*norm)
 #        m = min(5,m)

 #        # print "{:.5f}".format(m)

 #        if (am[-1] >= am[i_p_m]) :
 #            s_fc = 1
 #        else :
 #            s_fc = 1

 #        if prev_s != abs(s_fc)/s_fc :
 #            prev_s *= -1
 # #           f_fin = w/(2*pi)


def animate(t):

    global x, w, phi, ampl, f_fin, dt
    global s_fc
    global y_next, y_now, y_prev
    global cdtdx2
    global n
    global mic_t, mic
    global fft_spec, fft_axis, fft_axis_t, fft_peaks

    if chirp_type == "exp" :
        wnew = w * ( f_fin / f_ini )**(s_fc*dt/t_max)
    elif chirp_type == "lin" :
        # linear chirp gives trouble at low frequencies
        if ( w > pi ) | ( s_fc > 0 ) :
            wnew = w + s_fc*dt*2.0*pi*(f_fin-f_ini)/t_max;
        else :
            print "Low frequency bound reached."
            wnew = w
    else :
        print "ERROR: chirp_type not recognized."
        return

    phi = t*(w-wnew) + phi
    w = wnew

    boundary_y(t)
    step_y()

    if ( len(mic) > 10 ) :
        if ( n%30 == 0 ) :
            fft_norm = 2.0/len(mic_t)
            fft_spec = abs(fft_norm * fft.fft(mic))[:len(mic_t)/2]
            fft_axis = fft.fftfreq(len(mic_t), dt)[:len(mic_t)/2]
            # mph stands for minimum peak height
            fft_peaks = detect_peaks(fft_spec, mph=ampl/4.0)
    else :
        fft_spec = zeros(10)
        fft_axis = range(0,10)

    get_mic(mic_i, t)
    get_am(mic_i, t)

    n += 1

    if n%100 == 0 :
        print t

    # Change the value used to calculate the module according to how
    # offen do you want to refresh the plot.
    if n%1 == 0 :

        texts[0].set_text("t = " + "{:.4f}".format(t) + " s")
        texts[1].set_text("current emitted frequency = " + "{:.4f}".format(w/(2*pi)) + " Hz")
        texts[2].set_text("current step = " + str(n))
        # texts[3].set_text("s_fc = " + str(s_fc))
        peaks_str = ""
        for peak_i in fft_peaks :
            peaks_str += "{:.2f}, ".format(fft_axis[peak_i])
        peaks_str = "Peaks[Hz]: " + peaks_str[:-2]
        texts[9].set_text(peaks_str)

        # solution 'y_next' doesn't have the boundary points, that is
        # why 'y_now' is plotted.
        lines[0].set_ydata(y_now)
        # lines[1].set_ydata(pp(x,t))

        lines[2].set_xdata(mic_t)
        lines[2].set_ydata(mic)
        lines[3].set_xdata(mic_t)
        lines[3].set_ydata(am_min)
        lines[4].set_xdata(mic_t)
        lines[4].set_ydata(am_max)
        lines[5].set_xdata(mic_t)
        lines[5].set_ydata(am)

        scatters[1].set_offsets([[obs_i*dx],[y_now[obs_i]]])

        lines[6].set_xdata(f_ax)
        lines[6].set_ydata(am)

        lines[7].set_xdata(fft_axis)
        lines[7].set_ydata(5*fft_spec)

    y_prev = y_now.copy()
    y_now = y_next.copy()

    # You could return tuple(ax) but it makes the simulation slower.
    return tuple(lines) + tuple(texts) + tuple(scatters)

def init_ani() :

    global init_sim

    # When the window is resized init_ani is called again, but we
    # don't want our simulation to return to the initial state.
    if init_sim == 0 :
        init_y()

    init_sim += 1

    # It is neccesary to do this cleanup, otherwise objects get
    # overlaped.
    for line in lines :
        line.set_ydata(ma.array(line.get_xdata(), mask=True))

    for scatter in scatters :
        scatter.set_offsets([[],[]])

    for text in texts :
        text.set_text('')

    # costant values of the objects to be animated

    lines[6].set_linestyle('--')
    lines[7].set_linestyle('-')
    lines[7].set_marker('+')

    scatters[0].set_offsets([[mic_pos],[0]])

    texts[0].set_transform(ax[0].transAxes)
    texts[0].set_x(0.1)
    texts[0].set_y(0.9)
    texts[0].set_va('center')

    texts[1].set_transform(ax[0].transAxes)
    texts[1].set_x(0.4)
    texts[1].set_y(0.9)
    texts[1].set_va('center')

    texts[2].set_transform(ax[0].transAxes)
    texts[2].set_x(0.1)
    texts[2].set_y(0.75)
    texts[2].set_va('center')

    texts[3].set_transform(ax[0].transData)
    texts[3].set_x(mic_i*dx)
    texts[3].set_y(-2*ampl)
    texts[3].set_va('center')
    texts[3].set_text('measuring point, aka microphone')

    texts[4].set_transform(ax[0].transData)
    texts[4].set_x(obs_i*dx)
    texts[4].set_y(-2*ampl)
    texts[4].set_va('center')
    texts[4].set_text('obstacle')

    texts[5].set_transform(ax[1].transAxes)
    texts[5].set_x(0.1)
    texts[5].set_y(0.8)
    texts[5].set_va('center')
    texts[5].set_text('Mic signal, upper envelope, lower envelope, and amplitude.')

    texts[6].set_transform(ax[1].transAxes)
    texts[6].set_x(0.1)
    texts[6].set_y(0.65)
    texts[6].set_va('center')
    texts[6].set_text('Amplitude is the difference of the upper and lower ev.')

    texts[7].set_transform(ax[2].transAxes)
    texts[7].set_x(0.1)
    texts[7].set_y(0.8)
    texts[7].set_va('center')
    texts[7].set_text('Dashed: Amplitude as function of the current emitted freq.')

    texts[8].set_transform(ax[2].transAxes)
    texts[8].set_x(0.1)
    texts[8].set_y(0.65)
    texts[8].set_va('center')
    texts[8].set_text('Solid: 5 * FFT of the current signal')

    texts[9].set_transform(ax[2].transAxes)
    texts[9].set_x(0.1)
    texts[9].set_y(0.5)
    texts[9].set_va('center')

    # not used
    texts[10].set_transform(ax[2].transAxes)
    texts[10].set_x(0.1)
    texts[10].set_y(0.5)
    texts[10].set_va('center')

    return tuple(lines) + tuple(texts) + tuple(scatters)


# Create empty objects to be animated

lines = []
# first subplot
line = ax[0].plot(x,y_now)[0]   # this particular object cannot be empty :S
lines.append(line)
line = ax[0].plot([],[])[0]
lines.append(line)
# second subplot
line = ax[1].plot([],[])[0]
lines.append(line)
line = ax[1].plot([],[])[0]
lines.append(line)
line = ax[1].plot([],[])[0]
lines.append(line)
line = ax[1].plot([],[])[0]
lines.append(line)
# third subplot
line = ax[2].plot([],[])[0]
lines.append(line)
line = ax[2].plot([],[])[0]
lines.append(line)

scatters = []
scatter = ax[0].scatter([],[],s=10)
scatters.append(scatter)
scatter = ax[0].scatter([],[],s=10)
scatters.append(scatter)

texts = []
text = ax[0].text([], [], '')
texts.append(text)
text = ax[0].text([], [], '')
texts.append(text)
text = ax[0].text([], [], '')
texts.append(text)
text = ax[0].text([], [], '')
texts.append(text)
text = ax[0].text([], [], '')
texts.append(text)
text = ax[1].text([], [], '')
texts.append(text)
text = ax[1].text([], [], '')
texts.append(text)
text = ax[2].text([], [], '')
texts.append(text)
text = ax[2].text([], [], '')
texts.append(text)
text = ax[2].text([], [], '')
texts.append(text)
text = ax[2].text([], [], '')
texts.append(text)

ax[0].set_ylim( (-ampl*10,ampl*10) )
ax[0].set_xlim( 0,x[-1] )
ax[0].set_xticks( (0,L) )
ax[0].set_xlabel("x [m]",labelpad=-10)

ax[1].set_ylim( (-ampl*10,ampl*30) )
ax[1].set_xlim( 0,t_max )
ax[1].set_xlabel("time [s]",labelpad=-10)

ax[2].set_ylim( (0,ampl*30) )
ax[2].set_xlim( f_ini,f_fin )
ax[2].set_xlabel("frequency [Hz]",labelpad=0)

ani = animation.FuncAnimation(fig, animate, arange(t_min, t_max, dt),
                              repeat=False, init_func=init_ani,
                              interval=10, blit=True)

# mng = fig.canvas.manager
# # To maximize the animation window you can try the following
# # lines. Each of them could work, or not, depending on your
# # system. I didn't manage to make this work when generating a video
# # file.
# ##################
# # mng.frame.Maximize(True)
# # mng.window.showMaximized()
# # mng.window.state('zoomed')
# mng.resize(*mng.window.maxsize())  # tkagg
# ##################

if output == "fly" :
    show()
elif output == "file" :
    ani.save("string.mp4", writer="avconv", fps=25)
else :
    print "ERROR: output type not recognized."
