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

We consider a variational study of the square lattice antifer- romagnetic Heisenberg model with the Hamiltonian：

$$H=\sum_{\langle i,j\rangle}S_i^xS_j^x+S_i^yS_j^y+S_i^zS_j^z$$
![alt text](image-16.png)

---

![alt text](image-17.png)

---

![alt text](image-18.png)

---

![alt text](image-19.png)

---

![alt text](image-20.png)

---

![alt text](image-21.png)

---

![alt text](image-22.png)

---

![alt text](image-23.png)





























































