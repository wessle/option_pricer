from scipy.stats import norm
import numpy as np

# use the Black-Scholes formula to price a European call or put
def europ():
    CP = int(input("Enter 0 for call, 1 for put: "))
    S = float(input("Current price: "))
    K = float(input("Strike: "))
    tau = float(input("Expiration: "))
    sigma = float(input("Volatility: "))
    r = float(input("Risk-free interest rate: "))
    a = float(input("Continuous dividend rate: "))

    d1 = (1/(sigma*np.sqrt(tau)))*(np.log(S/K) + (r-a+(sigma**2/2))*tau)
    d2 = (1/(sigma*np.sqrt(tau)))*(np.log(S/K) + (r-a-(sigma**2/2))*tau)

    if CP == 0:
        price = S*np.exp(-a*tau)*norm.cdf(d1) - K*np.exp(-r*tau)*norm.cdf(d2)
    else:
        price = K*np.exp(-r*tau)*norm.cdf(-d2) - S*np.exp(-a*tau)*norm.cdf(-d1)

    return price
