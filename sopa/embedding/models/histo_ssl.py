# https://github.com/ozanciga/self-supervised-histopathology
import torch
from torchvision.models import resnet18


class HistoSSLFeatures(torch.nn.Module):
    def __init__(self):
        super().__init__()

        url = "https://github.com/ozanciga/self-supervised-histopathology/releases/download/tenpercent/tenpercent_resnet18.ckpt"
        # TODO: This checkpoint loading required pytorch_lightning it can be fixed by keeping only the weights and nothing extra from the checkpoint
        self.state = torch.hub.load_state_dict_from_url(url)
        self.model = resnet18(weights=None)
        self._load_model_weights()
        self.model.fc = torch.nn.Sequential()

    def _load_model_weights(self):
        state_dict = self.state["state_dict"]
        for key in list(state_dict.keys()):
            state_dict[key.replace("model.", "").replace("resnet.", "")] = state_dict.pop(key)

        model_dict = self.model.state_dict()
        state_dict = {k: v for k, v in state_dict.items() if k in model_dict}
        if state_dict == {}:
            print("No weight could be loaded..")
        model_dict.update(state_dict)
        self.model.load_state_dict(model_dict)

    def forward(self, x):
        return self.model(x)
