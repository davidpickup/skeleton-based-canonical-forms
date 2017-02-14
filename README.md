# Skeleton-Based Canonical Forms

This code is an implementation of the method presented in the following paper:

David Pickup, Xianfang Sun, Paul L. Rosin, Ralph R. Martin, "Skeleton-Based Canonical Forms for Non-Rigid 3D Shape Retrieval", Computational Visual Media, 2016.
(http://link.springer.com/article/10.1007/s41095-016-0045-5)


The implementation also contains code for computing the Laplace-Beltrami operator obtained from Mirela Ben-Chen, and redistributed with her permission, which is part of an implementation of the conformal factor measure presented in the following paper:

Ben-Chen, Mirela and Gotsman, Craig, Characterizing Shape Using Conformal Factors, Eurographics 2008 Workshop on 3D Object Retrieval, 2008, ISSN 1997-0463, 10.2312/3DOR/3DOR08/001-008.
(http://diglib.eg.org/handle/10.2312/3DOR.3DOR08.001-008)


The MDS implementation used is code released by Michael Bronstein.
(http://tosca.cs.technion.ac.il)


We request that if our code is used in your research that you cite the papers referenced above.

## Dependencies

The code requires the MDS implementation supplied by Michael Bronstein. The code can be downloaded
from either (http://tosca.cs.technion.ac.il/book/resources_sw.html), or is included as part of the
copy of this repository available at (https://sites.google.com/site/davidlpickup/software).

## Instructions

Download and extract the MDS code, as described above.

To compile the mex files and add the subdirectories to the Matlab path, please run the script "setup.m".

We have provided a sample mesh structure in the file "mesh.mat".

The function which computes our canonical form is:
"skeletonCanonicalForm.m"

To see a demo of our method, please run the script "demo.m".




THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.




