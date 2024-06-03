---
theme : "white"
transition: "fade"
center: false
height: 800
---
<style>
    .reveal h1, .reveal h2, .reveal h3, .reveal h4, .reveal h5 {
                  text-transform: none;
		  }
    .reveal p {
        text-align: left;
    }
    .reveal ul {
        display: block;
    }
    .reveal ol {
        display: block;
    }
    .reveal p:has(> img){
        text-align: center;
    }
    h3 {
        border-bottom: 2px solid yellow;
        padding: 10px;
    }
</style>


## Coner transfer matrix renormalization group algorithm

---

Matrix Product States (MPS) is a special form of tensor network, which is widely used to represent the wave function of 1D quantum many body systems， consider a 1D chain with N sites. We can see that:
![alt text](image-15.png)

---

For a 2D quantum many body system, we can use the Projected Entangled Pair States(peps) to represent the wave function, consider a 2D lattice with N sites. We can see that:

![alt text](image-14.png)

---

We consider an infinite projected entangled pair state (iPEPS) as the variational ansatz. The variational parameters are the elements in the iPEPS. Where s denotes the physical indices, and the remaining in- dices u, l, d, r are for virtual degrees of freedom of the bond dimension D. 

![alt text](image-13.png)

---

Consider 2 higher-order tensors $T^{s1,s2,s3,s4,s5,s6}$ and $W^{s1, s2, s3, s4, s5, s6}$, say that we want to compute the inner product of $T$ and $W$, if $W$=$T$this operation computes the norm of $T$. We want to compute:

![alt text](image-24.png)

---

We consider a variational study of the square lattice antifer- romagnetic Heisenberg model, and we want to find the ground state, minimize the expect energy: $\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$ 

We can see that the overlap of the iPEPS forms a tensor network, where the bulk tensor is the double layer tensor with bond dimension $d = D^2$

---

![alt text](image-16.png)

---

![alt text](image-23.png)

---

We can represent an infinite 2D square lattice with tensorsnetworks, consist of many many bulk tensors. We can see that:

![alt text](image-17.png)

---

We can define the edge tensors as a special combination of the bulk tensors, and we can see that:

![alt text](image-18.png)

---

We can also define the corner tensors as a special combination of the bulk tensors, and we can see that:

![alt text](image-19.png)

---

We want to solve this self-consistent equation iteratively.
![alt text](image-20.png)

---

![alt text](image-22.png)

---

![alt text](image-23.png)

---

![alt text](image-25.png)

---



























































