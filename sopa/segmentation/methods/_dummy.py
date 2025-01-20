import numpy as np


def dummy_method(**method_kwargs):
    """A method builder builder (i.e. it returns a segmentation function).
    Kwargs can be provided and used in the below function"""

    def segmentation_function(image: np.ndarray) -> np.ndarray:
        """A dummy example of a custom segmentation method
        that creates one cell (with a padding of 10 pixels).

        Args:
            image: An image of shape `(C, Y, X)`

        Returns:
            A mask of shape `(Y, X)` containing one cell
        """
        mask = np.zeros(image.shape[1:], dtype=int)

        # one cell, corresponding to value 1
        mask[10:-10, 10:-10] = 1  # squared shaped

        return mask

    return segmentation_function
