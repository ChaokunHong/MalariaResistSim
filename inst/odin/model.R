## Human Equations
deriv(S) <- -S * EIR_s * b_h * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) -
  S * EIR_R * b_h * (Phi * fT + Phi * (1 - fT) + (1 - Phi)) +
  Ts * rTs + As * rA + AR * rA + TR * rTR

deriv(Ds) <- S * EIR_s * b_h * Phi * (1 - fT) +
  EIR_s * b_h * As * Phi * (1 - fT) -
  Ds * rD

deriv(As) <- S * EIR_s * b_h * (1 - Phi) +
  Ds * rD -
  EIR_s * b_h * As * Phi * (1 - fT) -
  EIR_s * b_h * As * Phi * fT -
  As * rA

deriv(Ts) <- S * EIR_s * b_h * Phi * fT +
  EIR_s * b_h * As * Phi * fT -
  Ts * rTs

deriv(DR) <- S * EIR_R * b_h * Phi * (1 - fT) +
  EIR_R * b_h * AR * Phi * (1 - fT) -
  DR * rD

deriv(AR) <- S * EIR_R * b_h * (1 - Phi) +
  DR * rD -
  EIR_R * b_h * AR * Phi * (1 - fT) -
  EIR_R * b_h * AR * Phi * fT -
  AR * rA

deriv(TR) <- S * EIR_R * b_h * Phi * fT +
  EIR_R * b_h * AR * Phi * fT -
  TR * rTR

## Mosquito Equations
deriv(Sv) <- e - (Lambda_v_s + Lambda_v_r) * Sv - mu * Sv

delayed_Lambda_v_s_Sv <- delay(Lambda_v_s * Sv * exp(-mu * n), n)

deriv(Ev_s) <- Lambda_v_s * Sv - delayed_Lambda_v_s_Sv - mu * Ev_s
deriv(Iv_s) <- delayed_Lambda_v_s_Sv - mu * Iv_s

delayed_Lambda_v_r_Sv <- delay(Lambda_v_r * Sv * exp(-mu * n), n)

deriv(Ev_r) <- Lambda_v_r * Sv - delayed_Lambda_v_r_Sv - mu * Ev_r
deriv(Iv_r) <- delayed_Lambda_v_r_Sv - mu * Iv_r


## Initial conditions
initial(S) <- S0
initial(Ds) <- Ds0
initial(As) <- As0
initial(Ts) <- Ts0
initial(DR) <- DR0
initial(AR) <- AR0
initial(TR) <- TR0
initial(Sv) <- Sv0
initial(Ev_s) <- Ev_s0
initial(Iv_s) <- Iv_s0
initial(Ev_r) <- Ev_r0
initial(Iv_r) <- Iv_r0

## User-defined parameters
S0 <- user()
Ds0 <- user()
As0 <- user()
Ts0 <- user()
DR0 <- user()
AR0 <- user()
TR0 <- user()
EIR_s <- user()
EIR_R <- user()
Phi <- user()
fT <- user()
rD <- user()
rA <- user()
rTs <- user()
rTR <- user()
Sv0 <- user()
Ev_s0 <- user()
Iv_s0 <- user()
Ev_r0 <- user()
Iv_r0 <- user()
e <- user()
mu <- user()
n <- user()
Lambda_v_s <- user()
Lambda_v_r <- user()
b_h <- user()
