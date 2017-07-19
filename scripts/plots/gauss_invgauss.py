import sys, math

PI = 3.14159
scale = 0.915338301211 
shift = 18.4189782674

ig_sum = 0.0
g_sum = 0.0
mult_sum = 0.0

i = 0
while (True):
    try:
        lv_mean, lv_stdv, sd_mean, sd_stdv, sd_lambda = map(float, sys.stdin.readline().strip().split())
        ev_mean, ev_stdv = map(float, sys.stdin.readline().strip().split())
    except ValueError:
        break

    lv_mean = scale * lv_mean + shift
    #sd_mean = sd_mean

    g_prob = math.exp( -pow(ev_mean - lv_mean, 2) / (2 * pow(lv_stdv, 2)) ) / math.sqrt(2 * PI * pow(lv_stdv, 2))
    ig_prob = math.sqrt(sd_lambda / (2 * PI * pow(ev_stdv, 3)))  * math.exp(-sd_lambda * pow(ev_stdv - sd_mean, 2) / (2 * ev_stdv * pow(sd_mean, 2)));

    g_prob *= 1

    print ("%03d = %f\t%f\t%f" % (i, ig_prob*g_prob, g_prob, ig_prob))

    ig_sum += ig_prob
    g_sum += g_prob
    mult_sum += ig_prob*g_prob
    i += 1

print ("sum = %f\t%f\t%f" % (mult_sum, g_sum, ig_sum))
