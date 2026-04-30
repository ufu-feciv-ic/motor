# Uso de autovetores de flambagem como imperfeição inicial[cite: 1]
**Prof. Dr. Eduardo Vicente Wolf Trentini**

---

## Ideia central

Em uma análise não linear geométrica, uma estrutura perfeitamente reta e carregada de forma perfeitamente simétrica pode seguir uma trajetória idealizada. Para representar imperfeições geométricas iniciais, pode-se perturbar a geometria inicial usando uma forma modal de flambagem:

$$x_{0}^{imp}=x_{0}+e_{0}\phi_{n}$$

O vetor $\phi_{n}$ é uma versão normalizada do modo de flambagem, e $e_{0}$ é a amplitude física escolhida para a imperfeição.

---

## Sumário

| Seção | Tópico | Página |
| :--- | :--- | :--- |
| **1** | Objetivo | 3 |
| **2** | Por que não basta usar o autovetor de K? | 3 |
| **3** | Análise linear de referência | 3 |
| **4** | Matriz geométrica do elemento de pórtico 2D | 4 |
| **5** | Problema de flambagem linear | 5 |
| **6** | Redução para graus de liberdade livres | 5 |
| **7** | Conversão para um problema padrão | 6 |
| **8** | Normalização do modo | 7 |
| **9** | Alteração da geometria inicial | 7 |
| **10** | Algoritmo completo | 8 |
| **11** | Esboço de implementação com Eigen | 8 |
| **12** | Como escolher o modo | 9 |
| **13** | Cuidados numéricos e físicos | 10 |
| **13.1** | Apoios insuficientes | 10 |
| **13.2** | Sinal dos esforços normais | 10 |
| **13.3** | Escala da imperfeição | 10 |
| **13.4** | Não aplicar rotações como deslocamento | 10 |
| **13.5** | Magnitude da carga de referência | 10 |
| **14** | Exemplo conceitual: coluna comprimida | 11 |
| **15** | Resumo | 11 |

---

## 1. Objetivo

O objetivo deste texto é explicar como obter uma forma de imperfeição inicial para uma estrutura a partir de um problema de autovalor. A sequência usual é:

1. Montar a matriz de rigidez elástica $K$;
2. Aplicar uma carga de referência $F_{ref}$;
3. Resolver uma análise linear para obter $u_{ref}$;
4. Calcular os esforços normais de referência $N_{ref}$ nos elementos;
5. Montar a matriz de rigidez geométrica $K_{G}$;
6. Resolver o problema de flambagem linear;
7. Usar o primeiro modo de flambagem como imperfeição inicial.

A formulação apresentada é direcionada a estruturas reticuladas planas, especialmente elementos de pórtico 2D de Euler-Bernoulli. A mesma lógica geral aparece em outras formulações, mas a matriz geométrica do elemento muda.

---

## 2. Por que não basta usar o autovetor de K?

Se for resolvido o problema:

$$K\phi_{i}=\lambda_{i}\phi_{i}$$

Obtêm-se direções principais de rigidez da matriz elástica. O autovetor associado ao menor autovalor positivo pode indicar uma forma flexível da estrutura, mas não necessariamente representa o modo de instabilidade provocado pelo carregamento axial.

Na flambagem, o efeito importante é a redução da rigidez lateral causada por esforços normais de compressão. Esse efeito não está contido apenas em $K$. Ele entra pela matriz de rigidez geométrica $K_{G}$, montada a partir dos esforços normais existentes na estrutura. 

**Conclusão prática:** para imperfeição inicial associada à instabilidade, o mais adequado é usar o autovetor do problema de flambagem linear, e não simplesmente um autovetor de $K$.

---

## 3. Análise linear de referência

Primeiro, monta-se a rigidez elástica global da estrutura:

$$Ku=F$$

Escolhe-se uma carga de referência $F_{ref}$. Ela deve ter o mesmo padrão do carregamento que pode causar instabilidade. Por exemplo, se a estrutura será comprimida por cargas verticais nos nós superiores, essas cargas devem aparecer em $F_{ref}$. Pode ser inclusive o carregamento que se deseja aplicar à estrutura. A análise linear de referência é:

$$Ku_{ref}=F_{ref}$$

A partir do vetor $u_{ref}$, calcula-se o esforço normal de referência em cada elemento. Para um elemento de pórtico 2D, no sistema local, uma forma simples é:

$$N_{ref}^{(e)}=\frac{EA}{L}(u_{2,l}^{(e)}-u_{1,l}^{(e)})$$

Em que:
* $E$ é o módulo de elasticidade;
* $A$ é a área da seção;
* $L$ é o comprimento do elemento;
* $u_{1,l}^{(e)}$ e $u_{2,l}^{(e)}$ são os deslocamentos axiais locais dos nós do elemento.

Neste texto será adotada a convenção:

$$N>0 \text{ (tração)}, \quad N<0 \text{ (compressão)}$$

Com essa convenção, elementos comprimidos têm $N<0$ e reduzem a rigidez lateral da estrutura.

---

## 4. Matriz geométrica do elemento de pórtico 2D

Considere um elemento de pórtico plano de Euler-Bernoulli com graus de liberdade locais ordenados como:

$$d_{l}^{(e)}=\begin{bmatrix}u_{1}&v_{1}&\theta_{1}&u_{2}&v_{2}&\theta_{2}\end{matrix}^{T}$$

Uma forma usual da matriz geométrica local é:

$$k_{G,l}^{(e)} = \frac{N^{(e)}}{30L} \begin{bmatrix} 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & 36 & 3L & 0 & -36 & 3L \\ 0 & 3L & 4L^2 & 0 & -3L & -L^2 \\ 0 & 0 & 0 & 0 & 0 & 0 \\ 0 & -36 & -3L & 0 & 36 & -3L \\ 0 & 3L & -L^2 & 0 & -3L & 4L^2 \end{bmatrix}$$

Como a matriz acima está no sistema local do elemento, ela deve ser transformada para o sistema global:

$$k_{G}^{(e)}=T_{e}^{T}k_{G,l}^{(e)}T_{e}$$

Depois, cada matriz $k_{G}^{(e)}$ é assemblada na matriz global $K_{G}$, do mesmo modo que se faz com a matriz de rigidez elástica. 

**Sinal da matriz geométrica:** a matriz acima foi escrita com $N>0$ em tração e $N<0$ em compressão. Assim, para elementos comprimidos, contribui negativamente para a rigidez tangente lateral. Se outra convenção for usada, o sinal do problema de autovalor deve ser ajustado.

---

## 5. Problema de flambagem linear

A rigidez tangente linearizada em torno do estado de referência pode ser escrita como:

$$K_{T}(\lambda)=K+\lambda K_{G}$$

Em que $\lambda$ é um fator multiplicador da carga de referência. A perda de estabilidade linear ocorre quando existe um vetor não nulo $\phi$ tal que:

$$(K+\lambda K_{G})\phi=0$$

Essa equação pode ser reescrita como o problema de autovalor generalizado:

$$K\phi=-\lambda K_{G}\phi$$

O menor autovalor positivo $\lambda_{1}$ é o fator crítico de flambagem linear. Se $F_{ref}$ for a carga de referência, a carga crítica aproximada é:

$$F_{cr}\approx\lambda_{1}F_{ref}$$

O autovetor associado $\phi_{1}$ é o primeiro modo de flambagem.

---

## 6. Redução para graus de liberdade livres

Para resolver o problema de autovalor, é recomendável trabalhar apenas com os graus de liberdade livres. Separando os graus de liberdade em livres $f$ e restringidos $r$, escreve-se:

$$u=\begin{bmatrix}u_{f}\\ u_{r}\end{bmatrix}$$

Como os apoios impõem $u_{r}=0$, o problema de flambagem deve ser resolvido com as submatrizes livres:

$$K_{ff}\phi_{f}=-\lambda K_{Gff}\phi_{f}$$

Depois de calculado $\phi_{f}$, reconstrói-se o vetor global completo $\phi$, preenchendo com zero os graus de liberdade restringidos. 

> **Cuidado de implementação:** para análise estática linear, é comum impor apoios zerando linhas e colunas da matriz global. Para problema de autovalor, o procedimento mais limpo é montar diretamente as matrizes reduzidas $K_{ff}$ e $K_{Gff}$, contendo apenas os graus de liberdade livres.

---

## 7. Conversão para um problema padrão

Algumas bibliotecas resolvem diretamente o problema generalizado:

$$A\phi=\lambda B\phi$$

Nesse caso, pode-se usar:
* $A=K_{ff}$
* $B=-K_{Gff}$

Se a biblioteca disponível só resolver problemas padrão do tipo:

$$C\phi=\mu\phi$$

Então pode-se multiplicar a equação anterior por $K_{ff}^{-1}$:

$$K_{ff}^{-1}(-K_{Gff})\phi_{f}=\frac{1}{\lambda}\phi_{f}$$

Definindo:

$$C=K_{ff}^{-1}(-K_{Gff})$$

Tem-se:

$$C\phi_{f}=\mu\phi_{f} \quad \text{onde} \quad \mu=\frac{1}{\lambda}$$

Logo:

$$\lambda=\frac{1}{\mu}$$

Nesse formato, o menor $\lambda$ positivo corresponde ao maior $\mu$ positivo.

---

## 8. Normalização do modo

O autovetor não tem escala física. Se $\phi$ é autovetor, então $10\phi$, $0,01\phi$ ou $-\phi$ também são autovetores. Para usar o autovetor como imperfeição geométrica, deve-se escolher uma normalização. Para estruturas planas, uma opção simples é normalizar pelo maior deslocamento translacional:

$$\phi_{n}=\frac{\phi}{\max_{i}\sqrt{\phi_{x_{i}}^{2}+\phi_{y_{i}}^{2}}}$$

Assim, o maior deslocamento translacional nodal de $\phi_{n}$ passa a valer 1. A amplitude física da imperfeição pode ser definida por um valor $e_{0}$. Por exemplo:

$$e_{0}=\frac{L_{ref}}{1000}$$

Em que $L_{ref}$ pode ser o comprimento de uma barra, a altura de um pilar, o vão de uma estrutura, ou outra dimensão de referência adequada ao problema.

---

## 9. Alteração da geometria inicial

Depois da normalização, a geometria imperfeita é obtida por:

$$x_{i}^{imp}=x_{i}+e_{0}\phi_{x_{i},n}$$
$$y_{i}^{imp}=y_{i}+e_{0}\phi_{y_{i},n}$$

Em um pórtico plano, o autovetor costuma conter, para cada nó:

$$\begin{bmatrix}\phi_{x_{i}}&\phi_{y_{i}}&\phi_{\theta_{i}}\end{matrix}^{T}$$

Para modificar as coordenadas dos nós, usa-se normalmente apenas:
* $\phi_{x_{i}}$
* $\phi_{y_{i}}$

A componente rotacional $\phi_{\theta_{i}}$ pertence ao modo, mas não é uma coordenada geométrica nodal. Portanto, ela não precisa ser aplicada diretamente como alteração de geometria.

---

## 10. Algoritmo completo

A sequência de implementação pode ser resumida assim:

1. Montar $K$ elástica global.
2. Montar o vetor de carga de referência $F_{ref}$.
3. Identificar graus de liberdade livres.
4. Resolver $K_{ff}u_{ref,f}=F_{ref,f}$.
5. Reconstruir $u_{ref}$ global.
6. Para cada elemento - **Passo A**: transformar deslocamentos globais para locais;
7. Para cada elemento - **Passo B**: calcular $N_{ref}$ do elemento;
8. Para cada elemento - **Passo C**: montar $k_{G,l}^{(e)}$ local com $N_{ref}$;
9. Para cada elemento - **Passo D**: transformar $k_{G,l}^{(e)}$ para coordenadas globais;
10. Para cada elemento - **Passo E**: assemblar $k_{G}^{(e)}$ na matriz global $K_{G}$.
11. Reduzir $K$ e $K_{G}$ para os graus de liberdade livres.
12. Resolver $K_{ff}\phi=-\lambda K_{Gff}\phi$.
13. Escolher o menor $\lambda$ positivo, ou deixar o usuário escolher. Ou permitir a visualização dos modos e deixar o usuário escolher.
14. Reconstruir $\phi$ global.
15. Normalizar $\phi$ pelo maior deslocamento translacional.
16. Escolher $e_{0}$.
17. Criar geometria imperfeita: $x^{imp}=x+e_{0}\phi_{x}$; $y^{imp}=y+e_{0}\phi_{y}$.
18. Usar essa geometria como ponto inicial da análise não linear geométrica.

---

## 11. Esboço de implementação com Eigen

O trecho abaixo é apenas um esqueleto. Ele mostra a lógica, não uma implementação completa.
```cpp
// Matrizes globais completas
Eigen::MatrixXd K = AssembleElasticStiffness(model);
Eigen::VectorXd Fr = AssembleReferenceLoad(model);
Eigen::MatrixXd KG = Eigen::MatrixXd::Zero(K.rows(), K.cols());

// Lista de graus de liberdade livres
std::vector<int> freeDofs = BuildFreeDofList(model);

// Reducao para graus de liberdade livres
Eigen::MatrixXd Kff = ExtractSubmatrix(K, freeDofs, freeDofs);
Eigen::VectorXd Frf = ExtractSubvector(Fr, freeDofs);

// Analise linear de referencia
Eigen::VectorXd ur_f = Kff.ldlt().solve(Frf);
Eigen::VectorXd ur = ExpandToGlobalVector(ur_f, freeDofs, K.rows());

// Calculo dos esforcos normais e montagem de KG
for (Element& e : model.elements) {
    Eigen::VectorXd ue_global = ExtractElementDisplacements(ur, e);
    Eigen::VectorXd ue_local = e.T * ue_global;
    double N = e.E * e.A / e.L * (ue_local(3) - ue_local(0));
    Eigen::Matrix<double, 6, 6> kG_local = GeometricStiffness2DFrame(e.L, N);
    Eigen::Matrix<double, 6, 6> kG_global = e.T.transpose() * kG_local * e.T;
    AssembleElementMatrix(KG, kG_global, e.globalDofs);
}

Eigen::MatrixXd KGff = ExtractSubmatrix(KG, freeDofs, freeDofs);

// Problema padrao equivalente: C phi = mu phi
// C = Kff^{-1} * (-KGff), mu = 1/lambda
Eigen::MatrixXd C = Kff.ldlt().solve(-KGff);
Eigen::EigenSolver<Eigen::MatrixXd> solver(C);
// Selecionar maior mu positivo real; lambda = 1/mu.
// O autovetor associado e o modo de flambagem.