### Notation

### BI-2 Model

Definitions

$$
\begin{alignat}{4}
  & P_j & = & \lfloor p_j/\Delta \rfloor + 1, & \ & j \in J,\\
  & \pi_j & = & P_j - p_j/\Delta, & \ & j \in J,\\
  & D_j & = & \lfloor d_j/\Delta \rfloor + 1, & \ & j \in J,\\
  & \delta_j & = & D_j - d_j/\Delta , & \ & j \in J.
\end{alignat}
$$

Further definitions

$$
\begin{alignat}{4}
  & \pi_j^{\Delta} = \Delta \pi_j & = & \Delta P_j - p_j, & \ & j \in J,\\
  & \delta_j^{\Delta} = \Delta \delta_j & = & \Delta D_j - d_j, & \ & j \in J.
\end{alignat}
$$


#### Base formulation

$$
\begin{alignat}{2}
  \text{min}\ & \sum_{j \in J} \sum_{k \in K_j} \sum_{b \in [1, B-P_j+1]}
  c^z_{jb} z_{jbk} + c^u_{jb} u_{jbk},
  \\
  \text{s.t.}\ & \sum_{k \in K_j} \sum_{b \in [1, B-P_j+1]} z_{jbk} = 1, & \ &
  j \in J,
  \\
  & \sum_{j \in J} \sum_{k \in K_j} \sum_{a \in [b-P_j-k+2, b]} z_{jak}
  \leq 1, & \ & b \in [1, B],
  \\
  & \sum_{j \in J} \sum_{k \in K_j} (\Delta u_{jbk} - \Delta u_{j,b-P_j-k+1,k} &&\\
  & + (2\Delta-k\Delta-\pi_j^{\Delta}) z_{j,b-P_j-k+1,k} +
  \Delta \sum_{\substack{a \in [b-P_j-k+2, \\ b-1]}} z_{jak}) \leq \Delta, & \ & b
  \in [1, B],
  \\
  & ((1 - k)(\Delta - \pi_j^{\Delta}) + 1) z_{jbk} \leq \Delta u_{jbk}, & \
  & j \in J,\ k \in K_j,\ b \in [1, B-P_j+1],
  \\
  & \Delta u_{jbk} \leq (\Delta - k\pi_j^{\Delta}) z_{jbk}, & \ & j \in J,\ k \in K_j,\ b \in
  [1, B-P_j+1],
  \\
  & z_{jbk} \in \{0, 1\}, & \ & j \in J,\ k \in K_j,\ b \in [1, B-P_j+1].
\end{alignat}
$$

#### Capturing weighted tardiness

$$
\begin{alignat}{2}
  \text{min}\ & \sum_{j \in J} \sum_{k \in K} \sum_{b \in [D_j, B]}
  w_j \Delta T_{jbk},
  \\
  \text{s.t.}\ & \delta_j^{\Delta} z_{j D_j 0} - \Delta u_{j D_j 0} \leq \Delta T_{j
    D_j 0}, & \quad & j \in J \colon \Delta-\pi_j^{\Delta} < \delta_j^{\Delta},
  \\
  & \Delta T_{j D_j 0} \leq (\delta_j^{\Delta} - (\Delta - \pi_j^{\Delta})) z_{j D_j 0}, &
  \quad & j \in J \colon \Delta-\pi_j^{\Delta} < \delta_j^{\Delta},
  \\
  & \Delta T_{j D_j 1} = \delta_j^{\Delta} z_{j D_j 1} - \Delta u_{j D_j 1}, & \quad & j \in J \colon \Delta-\pi_j^{\Delta} < \delta_j^{\Delta},
  \\
  & \delta_j^{\Delta} z_{j D_j 1} - \Delta u_{j D_j 1} \leq \Delta T_{j D_j 1}, & \quad
  & j \in J \colon \delta_j^{\Delta} \leq \Delta-\pi_j^{\Delta},
  \\
  & \Delta T_{j D_j 1} \leq \delta_j^{\Delta} z_{j D_j 1}, & \quad & j \in J
  \colon \delta_j^{\Delta} \leq \Delta-\pi_j^{\Delta},
  \\
  & \Delta T_{jbk} = (\Delta b - \Delta D_j + \delta_j^{\Delta}) z_{jbk} - \Delta u_{jbk}, & \quad & j \in
  J,\ k \in K,\ _j \in [D_j + 1, B-P_j+1].
\end{alignat}
$$