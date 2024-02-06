import scipy.special

def get_legendre_coefficients(n):
    """
    Get the coefficients of the Legendre polynomial of degree n.

    Parameters:
    - n: The degree of the Legendre polynomial.

    Returns:
    - coeffs: An array of coefficients for the Legendre polynomial.
    """
    # Get the Legendre polynomial object of degree n
    legendre_poly = scipy.special.legendre(n)
    
    # The coefficients are stored in the .coef attribute
    coeffs = legendre_poly.coef
    return coeffs.tolist()

def get_Hermite_coefficients(n):
    """
    Get the coefficients of the Hermite polynomial of degree n.

    Parameters:
    - n: The degree of the Hermite polynomial.

    Returns:
    - coeffs: An array of coefficients for the Hermite polynomial.
    """
    # Get the Hermite polynomial object of degree n
    hermite_poly = scipy.special.hermitenorm(n)
    
    # The coefficients are stored in the .coef attribute
    coeffs = hermite_poly.coef
    return coeffs.tolist()

def get_legendre_quadrature_points(n):
    """
    Get the quadrature points of the Legendre polynomial of degree n.

    Parameters:
    - n: The degree of the Legendre polynomial.

    Returns:
    - points: A numpy array containing the quadrature points.
    """
    points, _ = scipy.special.roots_legendre(n)
    return points.tolist()

def get_Hermite_quadrature_points(n):
    """
    Get the quadrature points of the Hermite polynomial of degree n.

    Parameters:
    - n: The degree of the Hermite polynomial.

    Returns:
    - points: A numpy array containing the quadrature points.
    """
    points, _ = scipy.special.roots_hermite(n)
    return points.tolist()

def get_legendre_quadrature_weights(n):
    """
    Get the quadrature weights of the Legendre polynomial of degree n.

    Parameters:
    - n: The degree of the Legendre polynomial.

    Returns:
    - weights: A numpy array containing the quadrature weights.
    """
    _, weights = scipy.special.roots_legendre(n)
    return weights.tolist()

def get_Hermite_quadrature_weights(n):
    """
    Get the quadrature weights of the Hermite polynomial of degree n.

    Parameters:
    - n: The degree of the Hermite polynomial.

    Returns:
    - weights: A numpy array containing the quadrature weights.
    """
    _, weights = scipy.special.roots_hermite(n)
    return weights.tolist()
