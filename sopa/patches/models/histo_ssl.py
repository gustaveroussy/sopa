import torch


class HistoSSLFeatures(torch.nn.Module):
    def __init__(self):
        super().__init__()
        # weights retrieved from https://github.com/ozanciga/self-supervised-histopathology
        url = "https://nextcloud.centralesupelec.fr/s/Kn4zEGMdWgnj2ac/download/tenpercent_resnet12.pth"
        self.model = torch.hub.load_state_dict_from_url(url)
        self.model.fc = torch.nn.Sequential()
        self.output_dim = 512

    def forward(self, x):
        return self.model(x)
