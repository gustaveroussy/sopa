from spatialdata import SpatialData


class TranscriptsSegmentation:
    def __init__(self, sdata: SpatialData, method: callable):
        self.sdata = sdata
        self.method = method
