import numpy as np

def woodbury_update_covariance(sigma_te, c, lambda_gamma):
    """
    Computes the CE posterior covariance using the Woodbury rank-1 update.
    
    Formula:
        Sigma_CE = Sigma_TE - (lambda * Sigma_TE * c * c.T * Sigma_TE) / k
        where k = 1 + lambda * c.T * Sigma_TE * c
    
    Args:
        sigma_te (np.ndarray): The TE posterior covariance matrix (or diagonal).
        c (np.ndarray): The context vector (alm coefficients).
        lambda_gamma (float): The context coupling strength.
        
    Returns:
        np.ndarray: The CE covariance update (difference or full matrix).
    """
    # 1. Compute projection v = Sigma_TE * c
    # If sigma_te is 1D (diagonal approx), use element-wise multiply
    if sigma_te.ndim == 1:
        v = sigma_te * c
    else:
        v = sigma_te @ c
        
    # 2. Compute scalar denominator factor k
    # k = 1 + lambda * c.T * v
    c_dot_v = np.dot(c, v)
    k = 1.0 + lambda_gamma * c_dot_v
    
    # 3. Compute the update term
    # term = lambda/k * (v outer v)
    # This represents the "information loss" along the context direction.
    scale = lambda_gamma / k
    update = scale * np.outer(v, v)
    
    # In practice, we often return just the update or the diagonal of the new matrix
    # depending on memory constraints. Here we return the full rank-1 update matrix.
    return sigma_te - update

def calculate_mean_shift(mu_te, sigma_te, c, lambda_gamma, d_obs=None):
    """
    Computes the shift in the posterior mean: mu_CE - mu_TE.
    
    Formula:
        Delta_mu = lambda * (mu_TE.T * c) / (1 + lambda * c.T * Sigma_TE * c) * (Sigma_TE * c)
    
    Args:
        mu_te (np.ndarray): The TE posterior mean (Wiener filter).
        sigma_te (np.ndarray): The TE posterior covariance.
        c (np.ndarray): The context vector.
        lambda_gamma (float): Coupling strength.
        
    Returns:
        np.ndarray: The vector shift Delta_mu.
    """
    # 1. Compute projection v = Sigma_TE * c
    if sigma_te.ndim == 1:
        v = sigma_te * c
    else:
        v = sigma_te @ c
        
    # 2. Compute denominator k
    k = 1.0 + lambda_gamma * np.dot(c, v)
    
    # 3. Compute scalar overlap (mu_TE . c)
    # This measures how much the standard reconstruction looks like the context.
    overlap = np.dot(mu_te, c)
    
    # 4. Compute Shift
    # The shift is proportional to the overlap, directed along the "susceptibility" vector v.
    delta_mu = (lambda_gamma * overlap / k) * v
    
    return delta_mu
