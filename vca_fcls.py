import numpy as np
import matplotlib.pyplot as plt
from file_open import cuprite_data
from scipy.stats import norm
from pathlib import Path
from pysptools.abundance_maps import FCLS


def vca(Y, n_endmembers, snr_input=None):
    """
    Vertex Component Analysis (Nascimento & Dias, 2005).

    Parameters
    ----------
    Y : ndarray, shape (bands, pixels)
        Hyperspectral data matrix.
    n_endmembers : int
        Number of endmembers to extract.
    snr_input : float, optional
        SNR in dB. If None, it's estimated from the data.

    Returns
    -------
    E : ndarray, shape (bands, n_endmembers)
        Extracted endmember spectra.
    indices : ndarray, shape (n_endmembers,)
        Pixel indices in Y chosen as endmembers.
    """
    bands, pixels = Y.shape
    p = n_endmembers
    np.random.seed(42)

    # --- 1. Estimate SNR (decides projection strategy) ---
    y_mean = Y.mean(axis=1, keepdims=True)
    Y_centered = Y - y_mean

    Ud, _, _ = np.linalg.svd(Y_centered @ Y_centered.T / pixels)
    Ud = Ud[:, :p]

    if snr_input is None:
        proj = Ud.T @ Y_centered
        signal_power = np.mean(np.sum(Y ** 2, axis=0))
        noise_power = np.mean(np.sum(Y_centered ** 2, axis=0)) - np.mean(np.sum(proj ** 2, axis=0))
        snr_input = 10 * np.log10((signal_power - noise_power) / noise_power)

    snr_threshold = 15 + 10 * np.log10(p)

    # --- 2. Project data (two regimes depending on SNR) ---
    if snr_input > snr_threshold:
        # High SNR: project onto p principal components, keep affine structure
        x = Ud.T @ Y_centered
        u = x.mean(axis=1, keepdims=True)
        x_aug = x / (x * u).sum(axis=0, keepdims=True)  # normalize for numerical stability
    else:
        # Low SNR: project onto p-1 components and append a constant row
        Ud, _, _ = np.linalg.svd(Y @ Y.T / pixels)
        Ud = Ud[:, :p - 1]
        x = Ud.T @ Y
        c = np.sqrt(np.max(np.sum(x ** 2, axis=0)))
        x_aug = np.vstack([x, c * np.ones((1, pixels))])

    # --- 3. Iteratively find the p vertices (endmembers) ---
    indices = np.zeros(p, dtype=int)
    A = np.zeros((p, p))
    A[-1, 0] = 1  # seed the first direction

    for i in range(p):
        w = np.random.rand(p)
        # Project w onto the null space of the previously found vertices
        f = w - A @ np.linalg.pinv(A) @ w
        f = f / np.linalg.norm(f)

        # print(f.shape, x_aug.shape)
        proj_scores = f @ x_aug
        indices[i] = np.argmax(np.abs(proj_scores))
        A[:, i] = x_aug[:, indices[i]]

    E = Y[:, indices]
    return E, indices


def hfc(Y, false_alarm_rate=1e-4):
    """
    Harsanyi-Farrand-Chang virtual dimensionality estimator.
    Y: (bands, pixels)
    """
    bands, pixels = Y.shape

    # Correlation matrix (signal + noise)
    R = (Y @ Y.T) / pixels

    # Covariance matrix (mean-removed)
    Y_centered = Y - Y.mean(axis=1, keepdims=True)
    K = (Y_centered @ Y_centered.T) / pixels

    eig_R = np.sort(np.linalg.eigvalsh(R))[::-1]
    eig_K = np.sort(np.linalg.eigvalsh(K))[::-1]

    # Neyman-Pearson threshold based on desired false alarm rate
    tau = norm.ppf(1 - false_alarm_rate)

    p = 0
    for i in range(bands):
        sigma = np.sqrt(2 * (eig_R[i]**2 + eig_K[i]**2) / pixels)
        threshold = tau * sigma
        if (eig_R[i] - eig_K[i]) > threshold:
            p += 1

    return p


def hfc_details(Y, false_alarm_rate=1e-3, n_show=20):
    bands, pixels = Y.shape
    R = (Y @ Y.T) / pixels
    Y_centered = Y - Y.mean(axis=1, keepdims=True)
    K = (Y_centered @ Y_centered.T) / pixels

    eig_R = np.sort(np.linalg.eigvalsh(R))[::-1]
    eig_K = np.sort(np.linalg.eigvalsh(K))[::-1]

    tau = norm.ppf(1 - false_alarm_rate)

    for i in range(n_show):
        sigma = np.sqrt(2 * (eig_R[i]**2 + eig_K[i]**2) / pixels)
        threshold = tau * sigma
        diff = eig_R[i] - eig_K[i]
        print(f"{i:2d}  eig_R={eig_R[i]:.3e}  eig_K={eig_K[i]:.3e}  "
              f"diff={diff:.3e}  threshold={threshold:.3e}  "
              f"{'SIGNAL' if diff > threshold else 'noise'}")


cube_3d = cuprite_data()
# print(cube_3d.ndim)
cube_3d = cube_3d[:, :, 2:]
h, w, bands = cube_3d.shape
Y = cube_3d.reshape(-1, bands).T  # (bands, pixels)

# Center the data
Y_centered = Y - Y.mean(axis=1, keepdims=True)

# SVD on the centered data
U, S, Vt = np.linalg.svd(Y_centered, full_matrices=False)

# Plot the first ~30 singular values
n_show = 30
# plt.figure(figsize=(7, 5))
# plt.semilogy(range(1, n_show + 1), S[:n_show], 'o-')
# plt.xlabel('Component number')
# plt.ylabel('Singular value (log scale)')
# plt.title('SVD Singular Values (log scale)')
# plt.grid(alpha=0.3, which='both')
# plt.show()

Y_std = (Y) / Y.std(axis=1, keepdims=True)

hfc_details(Y_std)  # try on standardized data
p_estimate = hfc(Y_std)
print(f"HFC estimated number of endmembers: {p_estimate}")


p = p_estimate # your chosen endmember count, from hfc
E, idx = vca(Y, n_endmembers=p)
np.save('E.npy', E)

# Plot extracted endmember spectra
visualization_path = Path("visualizations")
visualization_path.mkdir(parents=True, exist_ok=True)
plot_path = visualization_path
graph_path = plot_path / "vca_extracted_endmembers.png"

if not graph_path.exists():
    for i in range(p):
        plt.plot(E[:, i], label=f'Endmember {i+1}') #clipped on edges to remove artifacts
    plt.xlabel('Band index')
    plt.ylabel('Reflectance')
    plt.legend()
    plt.title('VCA Extracted Endmembers')
    plt.savefig(plot_path / "vca_extracted_endmembers.png", dpi=300, bbox_inches="tight")


# VCA output E is (bands, p) — pysptools wants (p, bands)
E_pysptools = E.T  # shape (14, 188)

# cube_3d is already (rows, cols, bands) — exactly what pysptools expects
fcls = FCLS()
abundance_maps = fcls.map(cube_3d, E_pysptools)  # shape: (rows, cols, p)
np.save('abundance_maps.npy', abundance_maps)
print(abundance_maps.shape)  # should be (250, 190, 14)

fig, axes = plt.subplots(3, 5, figsize=(18, 10))
axes = axes.flatten()

for i in range(p_estimate):
    axes[i].imshow(abundance_maps[:, :, i], cmap='jet', vmin=0, vmax=1)
    axes[i].set_title(f'Endmember {i+1}')
    axes[i].axis('off')

# hide unused subplot
axes[14].axis('off')

plt.tight_layout()
plt.savefig(plot_path / "abundance_maps.png", bbox_inches="tight")
# plt.show()

total = abundance_maps.sum(axis=2)
print(total.min(), total.max(), total.mean())  # should all be very close to 1

Y = cube_3d.reshape(-1, cube_3d.shape[2]).T  # (bands, pixels)
A = abundance_maps.reshape(-1, abundance_maps.shape[2]).T  # (p, pixels)
Y_reconstructed = E @ A

rmse = np.sqrt(np.mean((Y - Y_reconstructed) ** 2))
print(Y.min(), Y.max(), Y.mean())
print(f"Reconstruction RMSE: {rmse:.4f}")