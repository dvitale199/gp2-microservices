from pydantic import BaseModel
from typing import Dict, List

class CarriersResponse(BaseModel):
    status: str
    message: str
    output_paths: Dict[str, str]