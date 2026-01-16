#include <LiquidCrystal.h>

#include <PriorityQueue.h>
#include <QueueList_Modified.h>

#include <DHT11.h>
#define ctsPin 2


//Variables

//Temperature Sensor
DHT11 dht11(4);
int temperature = 0;
int humidity = 0;
int result = 0;
int tempWaitTime = 5000;

//Touch Sensor
int touchPin = 2;
int screenVoltagePin = 9;
int ctsValue;
int ctsValue_past;
int touchWaitTime = 0;

// Screen
const int rs = 12, en = 13, d4 = 5, d5 = 10, d6 = 7, d7 = 8;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);

//Event loop + Priority Queue initialization
typedef struct Event {
  long executeTime;
  void (*run)(Event *event);
};

Event *event;
Event *touchEvent;
Event *tempEvent;
PriorityQueue<Event*> *queue;


long currentTime = 0;
long action = 0;
//Event functions

bool compareTime(Event *a, Event *b){
  return a->executeTime < b->executeTime;
}

//Check touch sensor
void checkTouch(Event *event){
  ctsValue = digitalRead(ctsPin);
  if (ctsValue == HIGH && ctsValue_past != ctsValue){
    digitalWrite(screenVoltagePin, !digitalRead(screenVoltagePin));
  }
  ctsValue_past = ctsValue;
  touchEvent = new Event();
  touchEvent->executeTime = millis();
  touchEvent->run = checkTouch;
  queue->push(touchEvent); 
  
}

//Check temperature
void checkTemperature(Event *event){
  temperature = 0;
  humidity = 0;
  result = dht11.readTemperatureHumidity(temperature, humidity);

  if (result == 0) {
      lcd.setCursor(0, 0);
      lcd.print(String("Temp: ") + temperature + String(" C"));

            
      Serial.print("Temperature: ");
      Serial.print(temperature);
      lcd.setCursor(0, 1);

      lcd.print(String("Humidity: ") + humidity + " %");


      Serial.print(" °C\tHumidity: ");
      Serial.print(humidity);
      Serial.println(" %");
  } else {
      Serial.println(DHT11::getErrorString(result));
  }
  
  tempEvent = new Event();
  tempEvent->executeTime = millis() + tempWaitTime;
  tempEvent->run = checkTemperature;
  queue->push(tempEvent); 

}




void setup() {
  Serial.begin(9600);
  lcd.begin(16, 2);
  pinMode(screenVoltagePin, OUTPUT);
  pinMode(ctsPin, INPUT);

  queue = new PriorityQueue<Event*>(compareTime);
  
  touchEvent = new Event();
  touchEvent->executeTime = 0;
  touchEvent->run = checkTouch;
  queue->push(touchEvent);

  tempEvent = new Event();
  tempEvent->executeTime = tempWaitTime;
  tempEvent->run = checkTemperature;
  queue->push(tempEvent); 
  

}


void loop() {
  currentTime = millis();
  if (!queue->isEmpty()){
    action = queue->peek()->executeTime;
    if (action <= currentTime){
      event = queue->pop();
      event->run(event);
      delete event;
    }
  }
}
