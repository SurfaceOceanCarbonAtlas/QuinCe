package uk.ac.exeter.QuinCe.web.datasets.plotPage;

import java.lang.reflect.Type;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonSerializationContext;
import com.google.gson.JsonSerializer;

public class PlotPageTableValueSerializer
  implements JsonSerializer<PlotPageTableValue> {

  @Override
  public JsonElement serialize(PlotPageTableValue src, Type typeOfSrc,
    JsonSerializationContext context) {

    JsonObject json = new JsonObject();

    json.addProperty("value", src.getValue());
    json.addProperty("used", src.getUsed());
    json.addProperty("qcFlag", src.getQcFlag().getFlagValue());
    json.addProperty("qcMessage", src.getQcMessage());
    json.addProperty("flagNeeded", src.getFlagNeeded());

    return json;
  }

}